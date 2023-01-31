import json
import re
from pathlib import Path

from Bio import Entrez
import idmapping

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
        self.gb_acc_id = self.ncbi_dir / "gb_acc_id"
        self.mapper_path = self.ncbi_dir / "mapper.json"

        with open(self.mapper_path, "r") as json_handle:
            self.mapper = json.load(json_handle)

        # let NCBI know who you are
        Entrez.email = "tobias.senoner@tum.de"
        Entrez.tool = "Biopython"
        Entrez.api_key = "cd66bc099e5133c5509d14b068ff8fa7bf08"

    def _get_json_path(self, acc_id: str) -> Path:
        return self.gb_acc_id / f"{acc_id}.json"

    def _search_json_file_locally(self, acc_id: str):
        json_file = None
        # if acc_id is in self.mapper
        if acc_id in self.mapper:
            json_file = self._get_json_path(acc_id=acc_id)
        # otherwise search number in crossreferences
        else:
            for gb_id, crossrefs in self.mapper.items():
                for _, db_uid in crossrefs.items():
                    if acc_id == db_uid:
                        json_file = self._get_json_path(acc_id=gb_id)
        return json_file

    def _save_mapper(self) -> None:
        with open(self.mapper_path, "w") as json_handle:
            json.dump(self.mapper, fp=json_handle, indent=4, sort_keys=True)

    def _save_entry(self, rec: dict) -> None:
        # save file
        acc_id = self.parse_acc_id(rec=rec)
        json_file = self._get_json_path(acc_id=acc_id)
        with open(json_file, "w") as json_handle:
            json.dump(rec, json_handle, indent=4)

        # update mapper
        crossrefs = self.parse_crossref_ids(rec=rec)
        self.mapper[acc_id] = crossrefs
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
        if db2use is None:
            raise Exception(
                f"No NCBI entries found for '{acc_id}' in dbs: {', '.join(dbs)}"
            )
        # make query to dedicated DB
        with Entrez.efetch(db=db2use, id=acc_id, **kwargs) as handle:
            for record in Entrez.parse(handle=handle):
                record = dict(record)
        return record

    @staticmethod
    def _load_json_file(json_file: Path) -> dict:
        with open(json_file, "r") as json_handle:
            data = json.load(json_handle)
        return data

    def get_record(self, gb_id: str) -> dict:
        json_file = self._search_json_file_locally(acc_id=gb_id)
        if json_file is None:
            data = self._fetch_entry(
                acc_id=gb_id,
                dbs=["nuccore", "protein"],
                rettype="gb",
                retmode="xml",
            )
            self._save_entry(rec=data)
        else:
            data = self._load_json_file(json_file=json_file)
        return data

    def get_crossref(self, gb_id: str) -> dict[str, str]:
        if gb_id in self.mapper:
            crossref = self.mapper[gb_id]
        else:
            rec = self.get_record(acc_id=gb_id)
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
        full_seq = None
        for feature in rec["GBSeq_feature-table"]:
            if feature["GBFeature_key"] == "CDS":
                for qual in feature["GBFeature_quals"]:
                    if qual["GBQualifier_name"] == "translation":
                        full_seq = qual["GBQualifier_value"]
            if not "GBFeature_quals" in feature:
                continue

        if (rec["GBSeq_moltype"] == "AA") and ("GBSeq_sequence" in rec):
            full_seq = rec["GBSeq_sequence"].upper()
        return full_seq

    def add_id2mapper(self, gb_id: str, new_id: dict[str, str]) -> None:
        self.mapper[gb_id].update(new_id)
        self._save_mapper()

    def get_uniprot_acc_id(self, rec: dict) -> tuple[str, str]:
        gb_id = self.parse_acc_id(rec=rec)
        crossfref = self.get_crossref(gb_id=gb_id)
        if "sp" in crossfref:
            acc_id, db = crossfref["sp"], "SP"
        elif "tr" in crossfref:
            acc_id, db = crossfref["tr"], "TR"
        elif "up_no" in crossfref:
            acc_id, db = None, None
        else:
            # Map IDs from GeneBank to UniProt using `idmapping`
            for from_db in ["EMBL-GenBank-DDBJ", "EMBL-GenBank-DDBJ_CDS"]:
                # get request: id_mapping
                job_id = idmapping.submit_id_mapping(
                    from_db=from_db, to_db="UniProtKB", ids=[gb_id]
                )
                if idmapping.check_id_mapping_results_ready(job_id):
                    link = idmapping.get_id_mapping_results_link(job_id)
                    results = idmapping.get_id_mapping_results_search(link)
                if "failedIds" in results:
                    continue
                result = results["results"][0]

                # parse DataBase
                entry_type = result["to"]["entryType"]
                if "TrEMBL" in entry_type:
                    db = "TR"
                elif "Swiss-Prot" in entry_type:
                    db = "SP"
                else:
                    raise Exception(f"Entry type isn't SP or TR, but {entry_type}")
                acc_id = result["to"]["primaryAccession"]
                acc_id, db = [acc_id, db]
                self.add_id2mapper(gb_id=gb_id, new_id={db.lower(): acc_id})
                break
            else:
                # if no entry was found
                acc_id = None
                self.add_id2mapper(gb_id=gb_id, new_id={"up_no": None})
        return acc_id, db