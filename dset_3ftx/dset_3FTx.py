# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on: Wed 18 Jan 2023 02:27:10
Description: Parse 3FTx dset and get uniprot ID by extraction or BLAST
Usage:       python 3FTx_dset.py

@author: tsenoner
"""
import argparse
import json
import re
from pathlib import Path

import idmapping
import ncbi_helper
import numpy as np
import pandas as pd
import uniprot_helper  # import BatchBlaster, get_tax_id
from Bio import Entrez
from pyfaidx import Fasta


def setup_arguments() -> argparse.Namespace:
    """Defines and parses required and optional arguments for the script"""
    parser = argparse.ArgumentParser(
        description="Merge & reshape PPIHP files output to a single CSV file"
    )

    parser.add_argument(
        "-ci",
        "--csv_file_in",
        required=True,
        type=str,
        help="Path to CSV file",
    )
    parser.add_argument(
        "-fi",
        "--fasta_file_in",
        required=True,
        type=str,
        help="Path to FASTA file",
    )
    parser.add_argument(
        "-co",
        "--csv_file_out",
        required=True,
        type=str,
        help="Path to CSV output file",
    )
    parser.add_argument(
        "-fo",
        "--fasta_file_out",
        required=True,
        type=str,
        help="Path to FASTA output file",
    )
    parser.add_argument(
        "-fg",
        "--genomic_fasta",
        required=True,
        type=str,
        help="FASTA file containing full sequences supported by genomic data",
    )
    parser.add_argument(
        "-fz",
        "--zhang_fasta",
        required=True,
        type=str,
        help="FASTA file with Bungarus multicinctus seqs from Zhang paper",
    )
    parser.add_argument(
        "-cr",
        "--ritu_csv",
        required=True,
        type=str,
        help="CSV file with mature sequences from Ritu paper",
    )
    parser.add_argument(
        "-t",
        "--taxon_mapper",
        required=True,
        type=str,
        help="Path to JSON file which maps species to taxon ids",
    )
    parser.add_argument(
        "-b",
        "--blast_dir",
        required=True,
        type=str,
        help="Path to directory where BLASTed JSON files are saved",
    )

    args = parser.parse_args()
    csv_in = Path(args.csv_file_in)
    fasta_in = Path(args.fasta_file_in)
    csv_out = Path(args.csv_file_out)
    fasta_out = Path(args.fasta_file_out)
    # additional fasta data
    genomic_fasta = Path(args.genomic_fasta)
    zhang_fasta = Path(args.zhang_fasta)
    ritu_csv = Path(args.ritu_csv)
    # paths to save tmp data
    taxon_mapper_file = Path(args.taxon_mapper)
    blast_dir = Path(args.out_dir)
    return (
        csv_in,
        fasta_in,
        csv_out,
        fasta_out,
        genomic_fasta,
        zhang_fasta,
        ritu_csv,
        taxon_mapper_file,
        blast_dir,
    )


def construct_df(csv_path, fasta_path):
    df = pd.read_csv(csv_path, sep=",")
    # remove rows with `snake genomic` as `Major group``
    # df = df.loc[~df["Major group"].isin(["snake genomic"]), :]

    # rename column
    df = df.rename(
        columns={
            "Embedding ID": "embedding_id",
            "Evolutionary order": "evolutionary_order",
            "Major group": "major_group",
            "Major taxon for the purposes of this study": "major_taxon",
            "Name in fasta": "fasta_id",
            "Name (existing or suggested)": "name",
            "Original fasta header": "original_fasta_header",
            "Preliminary cysteine group": "cystein_group",
            "Species": "species",
        }
    )

    # select columns to keep and add new ones
    cols2keep = [
        "cystein_group",
        "evolutionary_order",
        "fasta_id",
        "major_group",
        "major_taxon",
        "name",
        "original_fasta_header",
        "species",
    ]
    df = df[cols2keep]
    df["data_origin"] = "original"
    df["db"] = np.nan
    df["taxon_id"] = np.nan
    df["acc_id"] = np.nan

    # correct wrong taxa
    df.loc[df["species"] == "Micrurus_ tener", "species"] = "Micrurus tener"

    # add sequences to DataFrame from FASTA file
    seqs = (
        (header, str(seq).replace("-", ""))
        for header, seq in Fasta(str(fasta_path)).items()
    )
    df_seq = pd.DataFrame(seqs, columns=["fasta_id", "seq"])
    df = pd.merge(left=df, right=df_seq, on="fasta_id")

    # remove entry with same UniProt AccID but different/updated sequence
    df = df.drop(df[df["fasta_id"] == "Walterinnesia_aegyptia_C0HKZ8"].index)

    # add `genomic_id` + `data_origin` column
    df["genomic_id"] = df["original_fasta_header"].str.extract(
        r"^([a-zA-Z]{4}\w+) "
    )
    df.loc[df["genomic_id"].notna(), "data_origin"] = "genomic"
    print(f"- {len(df)} sequences from original analysis")
    return df


def add_genomic_full_seq(df: pd.DataFrame, fasta_path: Path) -> pd.DataFrame:
    genomic_fasta = Fasta(fasta_path)
    for header, seq in genomic_fasta.items():
        uid = header.split("_-_")[-1]
        # exclude sequences that do not start with M
        if not str(seq).startswith("M"):
            continue
        df.loc[df["genomic_id"] == uid, "full_seq"] = str(seq)

    nr_full_seq = len(df["full_seq"].dropna())
    print(
        f"- {nr_full_seq} full sequences information added by genomic supported"
        " alignment"
    )
    return df


def add_zhang_data(df: pd.DataFrame, fasta_path: Path) -> pd.DataFrame:
    """Add Bungarus multicinctus sequences from Zhang Zhi paper"""
    zhang_fasta = Fasta(fasta_path)

    data = {}
    for header, seq in zhang_fasta.items():
        seq = str(seq).replace("-", "")
        if header.startswith("Bmul"):
            data.setdefault("fasta_id", []).append(header)
            data.setdefault("full_seq", []).append(seq)
    df_zhang = pd.DataFrame(data=data)
    meta_data = {
        "species": "Bungarus multicinctus",
        "major_group": "3FTx",
        "major_taxon": "Elapidae",
        "data_origin": "paper_zhang",
    }
    df_zhang = df_zhang.assign(**meta_data)
    # df_zhang = df_zhang.drop_duplicates(subset="full_seq")
    df = pd.concat([df, df_zhang], ignore_index=True)
    print(
        f"- {len(df_zhang)} `Bungarus multicinctus` seqs added from Zhang paper"
    )
    return df


def add_ritu_data(
    df: pd.DataFrame,
    csv_path: Path,
    uniprot_collector: uniprot_helper.UniProtDataGatherer,
    ncbi_collector: ncbi_helper.NcbiDataGatherer,
) -> pd.DataFrame:
    """Add Drysdalia coronoides sequences from Ritu Chandna paper"""

    # --- load data
    df_ritu = pd.read_csv(csv_path)
    df_ritu = df_ritu.rename(
        columns={"uid": "fasta_id", "type": "cystein_group"}
    )
    df_ritu = df_ritu.assign(**{"data_origin": "paper_RituChandna"})

    # --- separate identifiers (gi_number and acc_id)
    df_ritu[["gi_number", "fasta_id"]] = df_ritu["fasta_id"].str.split(
        "|", n=1, expand=True
    )
    df_ritu["acc_id"] = df_ritu["gi_number"].str.extract(r"(^[^\d].+)")
    df_ritu["gi_number"] = df_ritu.loc[
        df_ritu["gi_number"].str.match(r"^\d"), "gi_number"
    ]

    # --- get metadata
    # get species names for UniProt entries
    for idx, row in df_ritu.loc[df_ritu["acc_id"].notna(), :].iterrows():
        data = uniprot_collector.get_entry(acc_id=row["acc_id"])
        species = data["organism"]["scientificName"]
        df_ritu.at[idx, "species"] = species
    # get NCBI entries using GI number and gather metadata
    for idx, row in df_ritu[df_ritu["gi_number"].notna()].iterrows():
        gi_nr = row["gi_number"]
        ncbi_rec = ncbi_collector.get_record(gb_id=gi_nr)
        acc_id, db = ncbi_collector.get_uniprot_acc_id(rec=ncbi_rec)
        species, _ = ncbi_collector.parse_taxon(rec=ncbi_rec)
        df_ritu.loc[idx, ["acc_id", "db", "species"]] = [acc_id, db, species]

    df = pd.concat([df, df_ritu], ignore_index=True)
    print(f"- {len(df_ritu)} `Drysdalia coronoides` seqs added from Ritu paper")
    return df


def create_taxon_mapper(taxas):
    taxon_mapper = {}
    for taxa in taxas:
        if taxa not in taxon_mapper:
            taxa_id = uniprot_helper.get_tax_id(taxa)
            if taxa_id is None:
                raise Exception(f"'{taxa}' not found")
            taxon_mapper[taxa] = taxa_id
    return taxon_mapper


def add_taxon_id(df: pd.DataFrame, taxon_mapper_file: Path) -> pd.DataFrame:
    # read in already mapped taxon
    if taxon_mapper_file.is_file():
        with open(taxon_mapper_file, "r") as json_file:
            taxon_mapper = json.load(json_file)
    else:
        taxon_mapper = dict()
    unknown_taxon_ids = df.loc[
        ~df["species"].isin(taxon_mapper.keys()), "species"
    ].unique()

    # get taxon_id that arn't in `taxon_mapper` yet
    novel_taxon_mapper = create_taxon_mapper(taxas=unknown_taxon_ids)
    taxon_mapper.update(novel_taxon_mapper)

    # update json file
    with open(taxon_mapper_file, "w") as json_file:
        json.dump(taxon_mapper, fp=json_file, indent=4)

    # add taxon ids to DataFrame
    df["taxon_id"] = df["species"].map(taxon_mapper)
    return df


def get_uniprot_acc_ids(df) -> tuple[list, dict]:
    # get UniProt ids
    for idx, row in df.iterrows():
        uid = "blank"
        header = row["original_fasta_header"]
        # 1. check for patterns that have no UniProt ID
        # - `3FTx_\d{3}`. E.g.: Cbivi_3FTx_000
        match1 = re.match(pattern=r".*(3FTx_\d{2,3})", string=header)
        # -`unigene\d*`. E.g.: Heterodon_nasicus_unigene14895
        match2 = re.match(pattern=r".*(unigene\d{4,6})", string=header)
        if match1 or match2:
            uid = None
        # 2. check for regular UniProt entry
        elif header.startswith("sp") or header.startswith("tr"):
            uid = header.split("|")[1]
        # 3. check for odd manually descriptive naming
        elif " " in header:
            uid = None
        # 4. check if last element of sequence is UID
        else:
            header_last_part = header.split("_")[-1]
            # UniProt regular expression for accession nummers: https://www.uniprot.org/help/accession_numbers
            uniprot_accession_pattern = r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
            match_accession = re.match(
                pattern=uniprot_accession_pattern, string=header_last_part
            )
            if match_accession:
                uid = header_last_part
            if not match_accession:
                uid = None
        df.loc[idx, "acc_id"] = uid
    print(f"- {df['acc_id'].count()} Acc IDs extracted ")
    print(f"- {df['acc_id'].isna().sum()} Acc IDs unknown")
    return df


def run_blast(
    df: pd.DataFrame,
    blast_dir: Path,
    uniprot_collector: uniprot_helper.UniProtDataGatherer,
) -> pd.DataFrame:
    # for entries having the full sequence use the full sequence to BLAST
    full_seq2blast = (
        df.loc[df["full_seq"].notna(), ["fasta_id", "full_seq", "taxon_id"]]
        .rename(columns={"full_seq": "seq"})
        .to_dict("records")
    )
    entries2blast = df.loc[
        df["full_seq"].isna(), ["fasta_id", "seq", "taxon_id"]
    ].to_dict("records")
    entries2blast.extend(full_seq2blast)

    # run BLASTp
    ncbi_blaster = uniprot_helper.UniProtBlaster(
        entries=entries2blast, out_dir=blast_dir
    )
    ncbi_blaster.run_batch()

    # extract ACC IDs
    for idx, row in df.iterrows():
        if not pd.isna(row["acc_id"]):
            continue

        # get accession IDs
        acc_id, db = ncbi_blaster.get_acc_id(
            fasta_id=row["fasta_id"], uniprot_collector=uniprot_collector
        )
        df.loc[idx, ["acc_id", "db"]] = [acc_id, db]

    return df


def get_uniprot_data(
    df: pd.DataFrame, uniprot_collector: uniprot_helper.UniProtDataGatherer
) -> pd.DataFrame:
    data = []
    for idx, row in df.iterrows():
        data = uniprot_collector.get_entry(acc_id=row["acc_id"])
        # check if a full sequence exists
        # - features: Signal + Chain
        # - different from existing sequence
        # - start with residue M
        features = set()
        if "features" in data:
            for feature in data["features"]:
                feature_type = feature["type"]
                features.add(feature_type)
        if set(["Signal", "Chain"]).issubset(features):
            json_seq = data["sequence"]["value"]
            if (json_seq != row["seq"]) and (json_seq.startswith("M")):
                row["full_seq"] = data["sequence"]["value"]

        # TODO: get the mature sequence add `Propeptide` if it is there
        # if "Propeptide" in features:
        #     print(acc_id, features)

        # get recomended name from entry
        for name_type, values in data["proteinDescription"].items():
            if name_type == "recommendedName":
                name = values["fullName"]["value"]
                df.loc[idx, "row"] = name

        # add other metadata
        metadata = [
            data["uniProtkbId"],
        ]
        df.loc[idx, ["acc_id", "entry", "prot_evi", "annot_score"]] = metadata


def manual_curation(df: pd.DataFrame) -> pd.DataFrame:
    # remove unneded column
    df = df.drop(columns="original_fasta_header")
    # remove duplicates
    len_before = len(df)
    df = df.drop_duplicates(subset=["acc_id", "species", "seq", "full_seq"])
    len_after = len(df)
    if len_before > len_after:
        print(f"{len_before-len_after} duplicate sequences were removed.")

    return df


def save_data(df: pd.DataFrame, csv_file: Path, fasta_file: Path) -> None:
    # move column `fasta_id` to the front
    col_data = df.pop("fasta_id")
    df.insert(loc=0, column="fasta_id", value=col_data)

    df.to_csv(csv_file, index=False)
    with open(fasta_file, "w") as fasta_handler:
        for _, row in df.iterrows():
            fasta_handler.write(f">{row['fasta_id']}\n")
            fasta_handler.write(f"{row['seq']}\n")


def main():
    (
        csv_in,
        fasta_in,
        csv_out,
        fasta_out,
        genomic_fasta,
        zhang_fasta,
        ritu_csv,
        taxon_mapper_file,
        blast_dir,
    ) = setup_arguments()

    df = construct_df(csv_path=csv_in, fasta_path=fasta_in)
    df = get_uniprot_acc_ids(df=df)
    df = add_genomic_full_seq(df=df, fasta_path=genomic_fasta)
    df = add_zhang_data(df=df, fasta_path=zhang_fasta)
    df = add_ritu_data(df=df, csv_path=ritu_csv)
    df = add_taxon_id(df=df, taxon_mapper_file=taxon_mapper_file)
    df = run_blast(df=df, blast_dir=blast_dir)
    df = manual_curation(df=df)
    save_data(df=df, csv_file=csv_out, fasta_file=fasta_out)


if __name__ == "__main__":
    main()
