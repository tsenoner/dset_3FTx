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

# --- patterns ---
# RefSeq ID: https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly
PATTERN_REFSEQ = r"([CMNPRXW]{2}_\d+)"  # {6,9})"
# GenBank ID, 3.4.6 ACCESSION Format: https://www.ncbi.nlm.nih.gov/genbank/release/current/
PATTERN_GENBANK = r"([A-Z]{1,4}\d{5,8})"
# UniProt ID confined: https://www.uniprot.org/help/accession_numbers
PATTERN_UNIPROT = (
    r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})"
)
CONFINEMENT1 = r"(?:^|[^A-Z0-9]|$)"  # everything but upper case char or number
CONFINEMENT2 = r"(?:^| |\.|_|$)"  # <start> <space> <dot> <underscore> <end>

# --- UniProt & NCBI API helpers ---
UNIPROT_DIR = Path("../data/uniprot_entries")
NCBI_DIR = Path("../data/ncbi_entries")
UNIPROT_COLLECTOR = uniprot_helper.UniProtDataGatherer(uniprot_dir=UNIPROT_DIR)
NCBI_COLLECTOR = ncbi_helper.NcbiDataGatherer(ncbi_dir=NCBI_DIR)


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


class OriginalDset:
    def __init__(
        self, csv_path: Path, fasta_path: Path, genomic_fasta_path: Path
    ) -> None:
        self.csv_path = csv_path
        self.fasta_path = fasta_path
        self.genomic_fasta_path = genomic_fasta_path

        df = self._basic_dset_preparation()
        df = self._extract_acc_ids(df=df)
        self.df = self._add_genomic_full_seq(df=df)

    def _basic_dset_preparation(self) -> pd.DataFrame:
        df = pd.read_csv(self.csv_path, sep=",")

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
        df[["uniprot_id", "refseq_id", "genbank_id"]] = np.nan

        # correct wrong taxa
        df.loc[df["species"] == "Micrurus_ tener", "species"] = "Micrurus tener"

        # add sequences to DataFrame from FASTA file
        seqs = (
            (header, str(seq).replace("-", ""))
            for header, seq in Fasta(str(self.fasta_path)).items()
        )
        df_seq = pd.DataFrame(seqs, columns=["fasta_id", "mature_seq"])
        df = pd.merge(left=df, right=df_seq, on="fasta_id")

        # remove entry with same UniProt AccID but different/old sequence
        # Walterinnesia_aegyptia_C0HKZ8 -> 353sp|C0HK
        df = df.drop(
            df[df["fasta_id"] == "Walterinnesia_aegyptia_C0HKZ8"].index
        )

        # add `genomic_id` + `data_origin` column
        df["genomic_id"] = df["original_fasta_header"].str.extract(
            r"^([a-zA-Z]{4}\w+) "
        )
        # df.loc[df["genomic_id"].notna(), "data_origin"] = "genomic"
        return df

    def _extract_acc_ids(self, df: pd.DataFrame) -> pd.DataFrame:
        # specific patterns to original DB
        pattern_3FTx = r"(3FTx_\d{2,3})$"
        pattern_unigene = r"(unigene\d{4,6})$"

        for idx, row in df.iterrows():
            header = row["original_fasta_header"]
            uniprot_match = re.search(
                CONFINEMENT1 + PATTERN_UNIPROT + CONFINEMENT1, header
            )
            refseq_match = re.search(PATTERN_REFSEQ, header)
            genbank_match = re.search(
                CONFINEMENT2 + PATTERN_GENBANK + CONFINEMENT2, header
            )
            FTx_match = re.search(pattern_3FTx, string=header)
            unigene_match = re.search(pattern_unigene, string=header)
            if uniprot_match:
                df.loc[idx, ["uniprot_id"]] = uniprot_match[1]
            elif refseq_match:
                df.loc[idx, ["refseq_id"]] = refseq_match[1]
            elif genbank_match:
                df.loc[idx, ["genbank_id"]] = genbank_match[1]
            elif FTx_match or unigene_match:
                continue
            else:
                if len(header.split(" ")) > 1:
                    # assemblies, could be searched but ignored for now
                    ncbi_search = header.split(" ", maxsplit=1)[1]

        print(
            f"{len(df)} entries: {df['uniprot_id'].count()} UniProt"
            f" IDs; {df['refseq_id'].count()} RefSeq IDs;"
            f" {df['genbank_id'].count()} GenBank IDs identified."
        )
        return df

    def _add_genomic_full_seq(self, df: pd.DataFrame) -> pd.DataFrame:
        genomic_fasta = Fasta(self.genomic_fasta_path)
        for header, seq in genomic_fasta.items():
            uid = header.split("_-_")[-1]
            # exclude sequences that do not start with M
            if not str(seq).startswith("M"):
                continue
            df.loc[df["genomic_id"] == uid, "full_seq"] = str(seq)

        nr_full_seq = len(df["full_seq"].dropna())
        print(
            f"- {nr_full_seq} full sequences information added by genomic"
            " supported alignment"
        )
        return df


class ZhangDset:
    def __init__(self, fasta_path: Path) -> None:
        self.fasta_path = fasta_path

        df = self._create_dset()
        self.df = self._extract_acc_ids(df=df)

    def _create_dset(self) -> pd.DataFrame:
        data = []
        for header, seq in Fasta(self.fasta_path).items():
            data.append([header, str(seq)])  # .replace("-", "")])
        df = pd.DataFrame(data=data, columns=["fasta_id", "mature_seq"])
        df["data_origin"] = "paper_zhang"

        # differentiate between full & mature sequences
        # full seq starts with M and in this alignment is before position 20
        full_seq_condition = (
            df["mature_seq"].str.replace("-", "").str.startswith("M")
        ) & (df["mature_seq"].str.find("M") <= 18)
        df["mature_seq"] = df["mature_seq"].str.replace("-", "")
        df.loc[full_seq_condition, "full_seq"] = df.loc[
            full_seq_condition, "mature_seq"
        ]
        df.loc[full_seq_condition, "mature_seq"] = None
        return df

    def _extract_acc_ids(self, df: pd.DataFrame) -> pd.DataFrame:
        # Elements found by GenBank pattern, but no entries when quering NCBI
        exclusion_lst = [
            "PDHV02000066",
            "PDHV02000188",
            "SOZL01001066",
            "SS00042983",
            "SS00017057",
        ]

        # extract Accession IDs
        for idx, row in df.iterrows():
            fasta_id = row["fasta_id"]
            uniprot_match = re.search(
                CONFINEMENT1 + PATTERN_UNIPROT + CONFINEMENT1, fasta_id
            )
            refseq_match = re.search(PATTERN_REFSEQ, fasta_id)
            genbank_match = re.search(PATTERN_GENBANK, fasta_id)
            if uniprot_match:
                df.loc[idx, ["uniprot_id"]] = uniprot_match[1]

            elif refseq_match:
                df.loc[idx, ["refseq_id"]] = refseq_match[1]

            elif (
                (genbank_match is not None)
                and (genbank_match[1] not in exclusion_lst)
                and ("scaffold" not in fasta_id)
            ):
                df.loc[idx, ["genbank_id"]] = genbank_match[1]

        print(
            f"{len(df)} entries: {df['uniprot_id'].count()} UniProt"
            f" IDs; {df['refseq_id'].count()} RefSeq IDs;"
            f" {df['genbank_id'].count()} GenBank IDs identified."
        )
        return df


class RituDset:
    def __init__(self, csv_path) -> None:
        self.csv_path = csv_path
        self.df = self._create_dset()

    def _create_dset(self) -> pd.DataFrame:
        df = pd.read_csv(self.csv_path)
        df = df.rename(
            columns={
                "uid": "fasta_id",
                "type": "ritu_class",
                "seq": "mature_seq",
            }
        )
        df = df.assign(**{"data_origin": "paper_ritu"})

        # --- separate identifiers (gi_number and acc_id)
        df[["gi_number", "fasta_id"]] = df["fasta_id"].str.split(
            "|", n=1, expand=True
        )
        df["uniprot_id"] = df["gi_number"].str.extract(r"(^[^\d].+)")
        df["gi_number"] = df.loc[df["gi_number"].str.match(r"^\d"), "gi_number"]
        print(
            f"{len(df)} entries: {df['uniprot_id'].count()} UniProt"
            f" IDs, {df['gi_number'].count()} GI numbers identified"
        )
        return df


class FrenchDset:
    def __init__(self, excel_path) -> None:
        self.excel_path = excel_path
        self.df = self._create_dset()

    def _create_dset(self) -> pd.DataFrame:
        # construct DataFrame + rename columns
        df = pd.read_excel(self.excel_path)
        df.drop(columns="SwissProt #", inplace=True)
        df = df.rename(
            columns={"Uniprot #": "uniprot_id", "Entry": "uniprot_entry"}
        )
        df.columns = [col.lower() for col in df.columns]
        cols2keep = [
            "uniprot_id",
            "uniprot_entry",
            "name",
            "receptor",
            "activity",
            "groups",
            "representative",
        ]
        df = df[cols2keep]

        # clean DataFrame up
        df.dropna(subset=["uniprot_id", "groups"], inplace=True)
        df.loc[df["representative"].isna(), "representative"] = False
        df.loc[df["representative"] == 1.0, "representative"] = True
        df["representative"] = df["representative"].astype(bool)
        df["groups"] = df["groups"].astype(int)
        df["data_origin"] = "french_guys"

        return df


def parse_uniprot_ids_file(uniprot_uids_files: list[Path]) -> pd.DataFrame:
    dfs = []
    for uniprot_uid_file in uniprot_uids_files:
        acc_ids = []
        if not uniprot_uid_file.is_file():
            raise Exception(f"{uniprot_uid_file} does not exist.")
        with open(uniprot_uid_file, "r") as handle:
            for line in handle:
                acc_id = line.strip()
                acc_ids.append(acc_id)
        df = pd.DataFrame(acc_ids, columns=["uniprot_id"])
        df["data_origin"] = uniprot_uid_file.stem
        dfs.append(df)
    df_uniprot = pd.concat(dfs)
    return df_uniprot


def map_ids2uniprot(df: pd.DataFrame) -> pd.DataFrame:
    # GI numbers to GenBank
    # df_map = NCBI_COLLECTOR.mapper
    # for idx, row in df.loc[df["gi_number"].notna(), :].iterrows():
    #     gb_id = df_map[df_map["gi"] == row["gi_number"]].values[0]
    #     df.loc[idx, "genbank_id"] = gb_id

    for idx, row in df.loc[df["gi_number"].notna(), :].iterrows():
        acc_id = ""
        NCBI_COLLECTOR.get_uniprot_acc_id(acc_id=acc_id)


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
    # TODO get family Rank of taxon
    # e.g. https://www.ebi.ac.uk/proteins/api/taxonomy/lineage/8633
    return df


def run_blast(
    df: pd.DataFrame,
    blast_dir: Path,
    uniprot_collector: uniprot_helper.UniProtDataGatherer,
) -> pd.DataFrame:
    # for entries having the full sequence use the full sequence to BLAST
    full_seq2blast = (
        df.loc[df["full_seq"].notna(), ["fasta_id", "full_seq", "taxon_id"]]
        .rename(columns={"full_seq": "mature_seq"})
        .to_dict("records")
    )
    entries2blast = df.loc[
        df["full_seq"].isna(), ["fasta_id", "mature_seq", "taxon_id"]
    ].to_dict("records")
    entries2blast.extend(full_seq2blast)

    # run BLASTp
    ncbi_blaster = uniprot_helper.UniProtBlaster(
        entries=entries2blast, out_dir=blast_dir
    )
    ncbi_blaster.run_batch()

    # extract ACC IDs
    for idx, row in df.iterrows():
        # if not pd.isna(row["acc_id"]):
        #     continue

        # get accession IDs
        acc_id = ncbi_blaster.get_acc_id(
            fasta_id=row["fasta_id"], uniprot_collector=uniprot_collector
        )
        df.loc[idx, "acc_id"] = acc_id

    return df


def get_uniprot_metadata(
    df: pd.DataFrame, uniprot_collector: uniprot_helper.UniProtDataGatherer
) -> pd.DataFrame:
    data = []
    for idx, row in df[df["acc_id"].notna()].iterrows():
        data = uniprot_collector.get_entry(acc_id=row["acc_id"])
        # get recomended name from entry
        for name_type, values in data["proteinDescription"].items():
            if name_type == "recommendedName":
                name = values["fullName"]["value"]
                df.loc[idx, "name"] = name

        # add other metadata
        metadata = [
            data["primaryAccession"],
            data["uniProtkbId"],
            uniprot_collector.parse_prot_evi(rec=data),
            data["annotationScore"],
        ]
        df.loc[idx, ["acc_id", "entry", "prot_evi", "annot_score"]] = metadata
    return df


def manual_curation(df: pd.DataFrame) -> pd.DataFrame:
    # TODO: remove sequences with a protein evidence level greater than 2
    # remove unneded column
    df = df.drop(columns="original_fasta_header")
    # TODO: merge duplicates to one record
    len_before = len(df)
    df = df.drop_duplicates(
        subset=["acc_id", "species", "mature_seq", "full_seq"]
    )
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

    # df = construct_df(csv_path=csv_in, fasta_path=fasta_in)
    # df = get_uniprot_acc_ids(df=df)
    # df = add_genomic_full_seq(df=df, fasta_path=genomic_fasta)
    # df = add_zhang_data(df=df, fasta_path=zhang_fasta)
    # df = add_ritu_data(df=df, csv_path=ritu_csv)
    # df = add_taxon_id(df=df, taxon_mapper_file=taxon_mapper_file)
    # df = run_blast(df=df, blast_dir=blast_dir)
    # df = manual_curation(df=df)
    # save_data(df=df, csv_file=csv_out, fasta_file=fasta_out)


if __name__ == "__main__":
    main()
