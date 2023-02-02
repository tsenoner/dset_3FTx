import re
from pathlib import Path

import ncbi_helper
import pandas as pd
import uniprot_helper
from pyfaidx import Fasta


def _get_metadata(
    ncbi_id: str,
    ncbi_collector: ncbi_helper.NcbiDataGatherer,
    uniprot_collector: uniprot_helper.UniProtDataGatherer,
):
    acc_id, db = ncbi_collector.get_uniprot_acc_id(gb_id=ncbi_id)
    if acc_id is not None:
        data = uniprot_collector.get_entry(acc_id=acc_id)
        species, taxon_id = uniprot_collector.parse_taxon(rec=data)
        full_seq, mature_peptide = uniprot_collector.parse_seq(rec=data)
    if (acc_id is None) or ((acc_id is not None) and (full_seq is None)):
        ncbi_rec = ncbi_collector.get_record(ncbi_id)
        species, taxon_id = ncbi_collector.parse_taxon(rec=ncbi_rec)
        full_seq = ncbi_collector.parse_seq(rec=ncbi_rec)
        mature_peptide = None
    return acc_id, db, species, taxon_id, full_seq, mature_peptide


def get_zhang_data(
    fasta_path: Path,
    uniprot_collector: uniprot_helper.UniProtDataGatherer,
    ncbi_collector: ncbi_helper.NcbiDataGatherer,
) -> pd.DataFrame:
    """Add Bungarus multicinctus sequences from Zhang Zhi paper"""
    zhang_fasta = Fasta(fasta_path)
    taxon_mapper = {
        "Bmul": "Bungarus multicinctus",
        "Cvir": "Crotalus viridis",
        "Dacu": "Deinagkistrodon acutus",
        "Hcur": "Hydrophis curtus",
        "Nnaj": "Naja naja",
        # "Pbiv": "Python bivittatus",
        # "Pgut": "Pantherophis guttatus",
        # "Tele": "Thamnophis elegans",
    }

    # --- patterns ---
    # UniProt ID confined: https://www.uniprot.org/help/accession_numbers
    pattern_uniprot = r"(?:_|^)([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(?:\.|_|$)"
    # RefSeq ID: https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly
    pattern_refseq = r"([CMNPRXW]{2}_\d{6,9})"
    # GenBank ID, 3.4.6 ACCESSION Format: https://www.ncbi.nlm.nih.gov/genbank/release/current/
    pattern_genbank = r"([A-Z]{1,4}\d{5,8})"
    # found by GenBank pattern, however there are no entries retrieved through API
    exclusion_lst = [
        "PDHV02000066",
        "PDHV02000188",
        "SOZL01001066",
        "SS00042983",
        "SS00017057",
    ]

    # --- map NCBI entries to UniProt
    ncbi_ids = []
    refseq_ids = []
    for header, _ in zhang_fasta.items():
        uniprot_match = re.search(pattern_uniprot, header)
        refseq_match = re.search(pattern_refseq, header)
        genbank_match = re.search(pattern_genbank, header)
        if uniprot_match:
            continue
        elif refseq_match:
            refseq_id = refseq_match[1]
            refseq_ids.append(refseq_id)
        elif (
            (genbank_match is not None)
            and (genbank_match[1] not in exclusion_lst)
            and ("scaffold" not in header)
        ):
            ncbi_id = genbank_match[0]
            ncbi_ids.append(ncbi_id)
    ncbi_collector.map_uniprot_acc_ids(
        ncbi_ids=refseq_ids, from_dbs=["RefSeq_Nucleotide", "RefSeq_Protein"]
    )
    ncbi_collector.map_uniprot_acc_ids(ncbi_ids=ncbi_ids)

    # --- create entries
    entries = []
    for header, seq in zhang_fasta.items():
        seq = str(seq).replace("-", "")
        uniprot_match = re.search(pattern_uniprot, header)
        refseq_match = re.search(pattern_refseq, header)
        genbank_match = re.search(pattern_genbank, header)
        # match in UniProt
        if uniprot_match:
            acc_id = uniprot_match[1]
            data = uniprot_collector.get_entry(acc_id=acc_id)
            species, taxon_id = uniprot_collector.parse_taxon(rec=data)
            full_seq, mature_seq = uniprot_collector.parse_seq(rec=data)
        # match in NCBI nuccore or protein DB
        elif refseq_match or (
            (genbank_match is not None)
            and (genbank_match[1] not in exclusion_lst)
            and ("scaffold" not in header)
        ):
            ncbi_id = (
                refseq_match[1] if genbank_match is None else genbank_match[1]
            )
            acc_id, db, species, taxon_id, full_seq, mature_seq = _get_metadata(
                ncbi_id=ncbi_id,
                ncbi_collector=ncbi_collector,
                uniprot_collector=uniprot_collector,
            )
        else:
            # BLASTp sequences
            taxa_abb = list(set(header.split("_")) & set(taxon_mapper.keys()))
            species = taxon_mapper[taxa_abb[0]]
            full_seq = seq
            acc_id, taxon_id, mature_seq = None, None, None

        entry = dict(
            fasta_id=header,
            acc_id=acc_id,
            db=db,
            full_seq=full_seq,
            seq=mature_seq,
            species=species,
            taxon_id=taxon_id,
            data_origin="paper_zhang",
        )
        entries.append(entry)

    df_zhang = pd.DataFrame(entries)
    df_zhang = df_zhang.dropna(subset=["mature_seq", "full_seq"], how="all")
    print(
        f"- {len(df_zhang)} `Bungarus multicinctus` seqs added from Zhang paper"
    )
    return df_zhang


def get_ritu_data(
    csv_path: Path,
    uniprot_collector: uniprot_helper.UniProtDataGatherer,
    ncbi_collector: ncbi_helper.NcbiDataGatherer,
) -> pd.DataFrame:
    """Return parsed Drysdalia coronoides sequences from Ritu Chandna paper"""

    # --- load data
    df_ritu = pd.read_csv(csv_path)
    df_ritu = df_ritu.rename(columns={"uid": "fasta_id", "type": "ritu_class"})
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
    ncbi_collector.map_uniprot_acc_ids(
        ncbi_ids=df_ritu[df_ritu["gi_number"].notna()]["gi_number"].to_list()
    )
    for idx, row in df_ritu[df_ritu["gi_number"].notna()].iterrows():
        gi_nr = row["gi_number"]
        acc_id, db, species, taxon_id, full_seq, mature_seq = _get_metadata(
            ncbi_id=gi_nr,
            ncbi_collector=ncbi_collector,
            uniprot_collector=uniprot_collector,
        )
        df_ritu.loc[
            idx,
            ["acc_id", "db", "species", "taxon_id", "full_seq", "mature_seq"],
        ] = [acc_id, db, species, taxon_id, full_seq, mature_seq]

    print(f"- {len(df_ritu)} `Drysdalia coronoides` seqs added from Ritu paper")
    return df_ritu

