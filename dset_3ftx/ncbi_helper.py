import json
import re
import sys
from pathlib import Path
from typing import Union
import pandas as pd

import idmapping
from Bio import Entrez

# documentation: https://www.nlm.nih.gov/dataguide/eutilities/utilities.html
# book: https://www.ncbi.nlm.nih.gov/books/NBK25501/

# RefSeq accession number format: https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly#
# UID for DBs: https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly
# retmode & rettype: https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly


# ----- Fetcher -----
# def entrez_efetch(
#     db: str, query: str, retmode: str = "xml", **kwargs
# ) -> list[dict]:
#     records = []
#     # with Entrez.efetch(db=db, id=identifier, rettype="asn.1",retmode="text", **kwargs) as handle:
#     #      records = handle.readlines()
#     with Entrez.efetch(db=db, id=query, retmode=retmode, **kwargs) as handle:
#         for rec in Entrez.parse(handle=handle):
#             records.append(dict(rec))
#     return records


# def get_ncbi_entry(ids: list[str], id_type: str, out_dir: Path) -> None:
#     ids = [uid for uid in ids if not (out_dir / f"{uid}.json").is_file()]
#     if ids:
#         for db_type in ["nucleotide", "protein"]:
#             for rec in entrez_efetch(
#                 db=db_type, query=",".join(ids), rettype="gb"
#             ):
#                 ids = parse_crossref_ids(rec=rec)
#                 uid = ids[id_type]
#                 ids.remove(uid)
#                 json_file = out_dir / f"{uid}.json"
#                 with open(json_file, "w") as json_handler:
#                     json.dump(rec, fp=json_handler, indent=4)


# def get_ncbi_data(db: str, uid: str, out_dir: Path) -> dict:
#     json_file = out_dir / f"{uid}.json"
#     if not json_file.is_file():
#         for data in entrez_efetch(db=db, query=uid, rettype="gb"):
#             with open(json_file, "w") as json_handle:
#                 json.dump(data, json_handle, indent=4)
#     else:
#         with open(json_file, "r") as json_handle:
#             data = json.load(json_handle)
#     return data


# ----- PARSERS -----
class NcbiDataGatherer:
    def __init__(self, ncbi_dir: Path) -> None:
        self.ncbi_dir = ncbi_dir
        self.gb_dir = self.ncbi_dir / "gb_entries"
        self.mapper_path = self.ncbi_dir / "mapper.csv"
        self.mapper = pd.read_csv(self.mapper_path)

        # let NCBI know who you are
        Entrez.email = "tobias.senoner@tum.de"
        Entrez.tool = "Biopython"
        Entrez.api_key = "cd66bc099e5133c5509d14b068ff8fa7bf08"

    def _save_mapper(self) -> None:
        self.mapper.to_csv(self.mapper_path, index=False)

    def _get_json_path(self, gb_id: str) -> Path:
        return self.gb_dir / f"{gb_id}.json"

    def _map_uid(self, query_id: str, to_id_type: str ="gb") -> str:
        to_id = None
        for col in self.mapper.columns:
            condition = (self.mapper[col] == query_id)
            if condition.sum():
                to_id = self.mapper.loc[condition, to_id_type].values[0]
                break
        return to_id

    def _search_json_file_locally(self, acc_id: str):
        json_file = None
        # if acc_id is in self.mapper
        if acc_id in self.mapper["gb"]:
            json_file = self._get_json_path(gb_id=gb_id)
        # otherwise search number in mapping file
        else:
            gb_id = self._map_uid(query_id=acc_id, to_id_type="gb")
            if gb_id is not None:
                json_file = self._get_json_path(gb_id=gb_id)
        return json_file

    def _save_entry(self, rec: dict) -> None:
        # save file
        acc_id = self.parse_acc_id(rec=rec)
        json_file = self._get_json_path(acc_id=acc_id)
        with open(json_file, "w") as json_handle:
            json.dump(rec, json_handle, indent=4)

        # update mapper
        crossrefs = self.parse_crossref_ids(rec=rec)
        self.mapper = pd.concat([self.mapper, pd.DataFrame([crossrefs])])
        self._save_mapper()

    def _fetch_entry(self, acc_id: str, dbs: list[str], **kwargs) -> list[dict]:
        # find out in which DB the entry can be found
        db2use = None
        with Entrez.egquery(term=acc_id) as handle:
            response = Entrez.read(handle=handle)
            for db_result in response["eGQueryResult"]:
                db_name, res_count = db_result["DbName"], db_result["Count"]
                if (res_count == "Error") or (res_count == "0"):
                    continue
                elif res_count == "1" and (db_name in dbs):
                    db2use = db_name
                # else:
                #     raise Exception(f"{db_name} found {res_count} entries.")
        if db2use is not None:
            # make query to dedicated DB
            with Entrez.efetch(db=db2use, id=acc_id, **kwargs) as handle:
                for record in Entrez.parse(handle=handle):
                    record = dict(record)
        else:
            print(
                f"No NCBI entries found for '{acc_id}' in dbs: {', '.join(dbs)}"
            )
            record = None
        return record

    @staticmethod
    def _load_json_file(json_file: Path) -> dict:
        with open(json_file, "r") as json_handle:
            data = json.load(json_handle)
        return data

    def get_record(self, gb_id: str) -> Union[dict, None]:
        json_file = self._search_json_file_locally(acc_id=gb_id)
        if json_file is None:
            data = self._fetch_entry(
                acc_id=gb_id,
                dbs=["nuccore", "protein"],
                rettype="gb",
                retmode="xml",
            )
            if data is not None:
                self._save_entry(rec=data)
        else:
            data = self._load_json_file(json_file=json_file)
        return data

    def get_crossref(self, gb_id: str) -> dict[str, str]:
        if gb_id in self.mapper["gb"]:
            crossref = self.mapper[self.mapper["gb"] == "gb"].to_dict("records")[0]
        else:
            rec = self.get_record(gb_id=gb_id)
            crossref = self.parse_crossref_ids(rec=rec)
        return crossref

    @staticmethod
    def parse_crossref_ids(rec: dict) -> dict[str, str]:
        """Extract crossreferencing IDs and parse them into a dict"""
        seq_ids = dict()
        for seq_id in rec["GBSeq_other-seqids"]:
            match = re.match(
                r"^([a-z]+)\|([\w:]+)(?<=\.\d\|)?", seq_id
            ).groups()
            seq_ids[match[0]] = match[1]
        return seq_ids

    @staticmethod
    def parse_acc_id(rec: dict) -> str:
        acc_id = rec["GBSeq_primary-accession"]
        return acc_id

    @staticmethod
    def _loop_quals(rec: dict):
        for feature in rec["GBSeq_feature-table"]:
            if not "GBFeature_quals" in feature:
                continue
            for qual in feature["GBFeature_quals"]:
                yield qual

    def parse_taxon(self, rec: dict) -> tuple[str, str]:
        species, taxon_id = None, None
        for qual in self._loop_quals(rec=rec):
            if qual["GBQualifier_name"] == "organism":
                species = qual["GBQualifier_value"]
            if qual["GBQualifier_name"] == "db_xref":
                taxon_id = qual["GBQualifier_value"].split(":")[1]
        return species, taxon_id

    def parse_seq(self, rec: dict) -> str:
        # TODO: quality check when sequence can be trusted
        entry_type = rec["GBSeq_moltype"]
        # if entry_type not in ["mRNA", "AA"]:
        #     print(entry_type, rec["GBSeq_length"], rec["GBSeq_primary-accession"])
        full_seq = None
        prot_id = None
        for feature in rec["GBSeq_feature-table"]:
            if feature["GBFeature_key"] == "CDS":
                for qual in feature["GBFeature_quals"]:
                    if qual["GBQualifier_name"] == "translation":
                        full_seq = qual["GBQualifier_value"]
                    elif qual["GBQualifier_name"] == "protein_id":
                        prot_id = qual["GBQualifier_value"]
            if not "GBFeature_quals" in feature:
                continue

        if (rec["GBSeq_moltype"] == "AA") and ("GBSeq_sequence" in rec):
            full_seq = rec["GBSeq_sequence"].upper()

        if (full_seq is not None) and (entry_type not in ["mRNA", "AA"]) and (prot_id is None):
            print(
                entry_type, rec["GBSeq_length"], rec["GBSeq_primary-accession"], prot_id
            )
        elif (full_seq is None) and (entry_type in ["mRNA", "AA"]):
            pass
        return full_seq

    def _update_mapper(self, gb_id: str, new_id: dict[str, str]) -> None:
        self.mapper[gb_id].update(new_id)
        self._save_mapper()

    def map_uniprot_acc_ids(self, ncbi_ids: list[str], from_dbs: list[str] = ["EMBL-GenBank-DDBJ", "EMBL-GenBank-DDBJ_CDS"]) -> None:
        # check which entries don't have a UniProt mapping (sp, tr, no_uniprot)
        gb_ids = list()
        for ncbi_id in ncbi_ids:
            ncbi_rec = self.get_record(ncbi_id)
            gb_id = self.parse_acc_id(rec=ncbi_rec)
            crossref = self.get_crossref(gb_id=gb_id)
            if not (set(["sp", "tr", "no_uniprot"]) & set(crossref.keys())):
                gb_ids.append(gb_id)

        # Map IDs from GeneBank to UniProt using `idmapping`
        if gb_ids:
            for from_db in from_dbs:
                # get request: id_mapping
                job_id = idmapping.submit_id_mapping(
                    from_db=from_db, to_db="UniProtKB", ids=gb_ids
                )
                if idmapping.check_id_mapping_results_ready(job_id):
                    link = idmapping.get_id_mapping_results_link(job_id)
                    results = idmapping.get_id_mapping_results_search(link)
                for result in results["results"]:
                    entry_type = result["to"]["entryType"]
                    if "TrEMBL" in entry_type:
                        db = "TR"
                    elif "Swiss-Prot" in entry_type:
                        db = "SP"
                    else:
                        raise Exception(f"Not SP or TR, but {entry_type}")
                    acc_id = result["to"]["primaryAccession"]
                    gb_id = result["from"]
                    if "no_uniprot" in self.mapper[gb_id]:
                        self.mapper[gb_id].pop("no_uniprot")
                    self._update_mapper(gb_id=gb_id, new_id={db.lower(): acc_id})
                gb_ids = results["failedIds"] if "failedIds" in results else []
            for failed_entry in gb_ids:
                # if no UniProt entry was found
                self._update_mapper(gb_id=failed_entry, new_id={"no_uniprot": None})

    def get_uniprot_acc_id(self, acc_id: str) -> tuple[str, str]:
        gb_id = self._map_uid(query_id=acc_id, to_id_type="gb")
        crossref = self.get_crossref(gb_id=gb_id)
        if "sp" in crossref:
            acc_id, db = crossref["sp"], "SP"
        elif "tr" in crossref:
            acc_id, db = crossref["tr"], "TR"
        elif "no_uniprot" in crossref:
            acc_id, db = None, None
        else:
            print(crossref)
            raise Exception(f"Map {gb_id} first to UniProt.")

        return acc_id, db
