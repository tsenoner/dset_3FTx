import json
import sys
from pathlib import Path
from typing import Union

import requests

# Documentation: https://www.uniprot.org/help/api
UNIPROT_API = "https://rest.uniprot.org"

# Documentation: https://www.ebi.ac.uk/proteins/api/doc/
PROTEINS_API = "https://www.ebi.ac.uk/proteins/api"

# https://www.ebi.ac.uk/proteins/api/doc/
TAXA_API = "https://www.ebi.ac.uk/proteins/api/taxonomy"

# Helper function to download data
def get_url(url, **kwargs):
    response = requests.get(url, **kwargs)

    if not response.ok:
        print(response.text)
        response.raise_for_status()
        sys.exit()

    return response


def get_taxa_id(taxa: str) -> Union[None, str]:
    """For a given taxa name, search and return its taxa ID.

    Return the taxa ID as an integer if it exists, else return None.
    """
    url = f"{TAXA_API}/name/{taxa}"
    response = get_url(url)
    try:
        taxonomies = response.json()["taxonomies"]
    except requests.HTTPError as err:
        txt = err
        if str(txt).startswith("404 Client Error: Not Found for url"):
            taxonomies = None
        else:
            raise Exception(err)

    tax_id = None
    if taxonomies is None:
        tax_id = None
    elif len(taxonomies) == 1:
        tax_id = taxonomies[0]["taxonomyId"]
    else:
        raise Exception(f"{taxa} found more than one taxa.")
    return tax_id


def get_taxa_rank(taxa_id: Union[str, int], rank: str = "family") -> str:
    url = f"{TAXA_API}/lineage/{taxa_id}"
    response = requests.get(url=url)
    rank_name = "None"
    for taxonomy in response.json()["taxonomies"]:
        if taxonomy["rank"] == rank:
            rank_name = taxonomy["scientificName"]
            break
    return rank_name


def load_data(json_path: Path) -> dict:
    with open(json_path, "r") as json_handle:
        data = json.load(json_handle)
    return data


class UniProtDataGatherer:
    def __init__(self, uniprot_dir: Path) -> None:
        self.uniprot_dir = uniprot_dir

    def _get_json_path(self, acc_id: str) -> Path:
        return self.uniprot_dir / f"{acc_id}.json"

    def _save_entry(self, rec: dict) -> None:
        acc_id = rec["primaryAccession"]
        json_file = self._get_json_path(acc_id=acc_id)
        with open(json_file, "w") as json_handle:
            json.dump(rec, fp=json_handle, indent=4)

    def _search_acc_id_locally(self, acc_id: str) -> Union[dict, None]:
        json_data = None
        json_file = self._get_json_path(acc_id=acc_id)
        if json_file.is_file():
            json_data = load_data(json_path=json_file)
        return json_data

    def get_entry(self, acc_id: str, format_type: str = "json") -> dict:
        """Get the result of a single UniProt entry by accession ID"""
        json_data = self._search_acc_id_locally(acc_id=acc_id)
        if json_data is None:
            url = f"{UNIPROT_API}/uniprotkb/{acc_id}.{format_type}"
            json_data = get_url(url=url).json()
            self._save_entry(rec=json_data)
        return json_data

    def parse_taxon(self, rec: dict) -> tuple[str, str]:
        """Parses species scientific name and taxon id"""
        species = rec["organism"]["scientificName"]
        taxon_id = rec["organism"]["taxonId"]
        return species, taxon_id

    @staticmethod
    def _get_start_end_idx(seq_idx):
        """Sort indices and remove any consecutive idx that increase by one

        Those indices mean that the sequence is consecutive
        """
        seq_idx = sorted(seq_idx) + [None]
        new_idx = []
        skip_flag = False
        for idx, jdx in zip(seq_idx, seq_idx[1:]):
            if (idx + 1 == jdx) or (idx == jdx):
                skip_flag = True
            elif skip_flag:
                skip_flag = False
            else:
                new_idx.append(idx)
        return new_idx

    def parse_name(self, rec: dict) -> Union[str, None]:
        # get recomended name from entry
        name = None
        prot_descr = rec["proteinDescription"]
        if "recommendedName" in prot_descr:
            name = prot_descr["recommendedName"]["fullName"]["value"]
        # elif "submissionNames" in prot_descr:
        #     name = prot_descr["submissionNames"][0]["fullName"]["value"]
        return name

    def parse_seq(self, rec: dict) -> tuple[Union[str, None], Union[str, None]]:
        """return full_seq (Signal + mature_seq) and mature_seq (Chain + Propeptide)

        full_seq = Signal + mature_seq
        mature_seq = (Propeptide) + Chain + (Propeptide)
        """
        # UniProt options for 3FTx
        # signal + chain: https://www.uniprot.org/uniprotkb/P19959/entry#ptm_processing
        #     propeptide: https://www.uniprot.org/uniprotkb/A0S865/entry#ptm_processing
        #                 https://www.uniprot.org/uniprotkb/Q8IV16/entry#ptm_processing
        #          chain: https://www.uniprot.org/uniprotkb/A0A4P1LYC9/entry#ptm_processing
        #        peptide: https://www.uniprot.org/uniprotkb/C0HJT4/entry#ptm_processing
        #           None: https://www.uniprot.org/uniprotkb/A0A6P9C7G6/entry#ptm_processing
        #   non-terminal: https://www.uniprot.org/uniprotkb/A1IVR8/entry#sequences

        # --- get sequence feature indices
        non_term_pos, non_adj_pos = None, None
        seq_featurs = {}
        seq_value = rec["sequence"]["value"]
        if "features" in rec:
            for feature in rec["features"]:
                feature_type = feature["type"]
                if feature_type in ["Chain", "Peptide", "Propeptide", "Signal"]:
                    start = feature["location"]["start"]["value"]
                    end = feature["location"]["end"]["value"]
                    seq_featurs.setdefault(feature_type, []).extend(
                        [start, end]
                    )
                elif feature_type == "Non-terminal residue":
                    non_term_pos = feature["location"]["start"]["value"]
                elif feature_type == "Non-adjacent residue":
                    non_adj_pos = feature["location"]["start"]["value"]

        acc_id = rec["primaryAccession"]
        # if acc_id == "A1IVR8":
        #     print(seq_featurs)
        #     print(non_term_pos, non_adj_pos)

        # --- get mature peptide
        if "Propeptide" in seq_featurs:
            mature_pep_idx = seq_featurs["Propeptide"] + seq_featurs["Chain"]
            mature_pep_idx = self._get_start_end_idx(seq_idx=mature_pep_idx)
        elif "Chain" in seq_featurs:
            mature_pep_idx = seq_featurs["Chain"]
        elif "Peptide" in seq_featurs:
            mature_pep_idx = seq_featurs["Peptide"]
        elif seq_featurs == {}:
            mature_pep_idx = None
        else:
            print(f"Sequence features: {seq_featurs}")
            raise Exception(
                f"Odd seq features for {rec['primaryAccession']}"
            )  #

        if (
            (mature_pep_idx is None)
            or (
                (non_term_pos is not None)
                and (mature_pep_idx[0] <= non_term_pos >= mature_pep_idx[1])
            )
            or (
                (non_adj_pos is not None)
                and (non_adj_pos[0] <= non_adj_pos >= non_adj_pos[1])
            )
        ):
            mature_peptide = None
        elif len(mature_pep_idx) == 2:
            mature_peptide = seq_value[
                mature_pep_idx[0] - 1 : mature_pep_idx[1]
            ]
        else:
            print(f"Sequence features: {seq_featurs}")
            raise Exception(f"Odd seq features for {rec['primaryAccession']}")

        # --- get full sequence
        if (
            (mature_pep_idx is None)
            or (non_term_pos is not None)
            or (non_adj_pos is not None)
            or ("Signal" not in seq_featurs)
        ):
            full_seq = None
        else:
            full_seq_idx = seq_featurs["Signal"] + mature_pep_idx
            full_seq_idx = self._get_start_end_idx(seq_idx=full_seq_idx)
            if len(full_seq_idx) == 2:
                full_seq = seq_value[full_seq_idx[0] - 1 : full_seq_idx[1]]
                if not full_seq.startswith("M"):
                    full_seq = None
            else:
                print(f"Sequence features: {seq_featurs}")
                print(f"Sequence indices: {full_seq_idx}")
                raise Exception(
                    f"Odd seq features for {rec['primaryAccession']}"
                )

        return full_seq, mature_peptide

    def parse_prot_evi(self, rec: dict) -> int:
        """Return the protein evidence level in UniProt

        1. Experimental evidence at protein level
        2. Experimental evidence at transcript level
        3. Protein inferred from homology
        4. Protein predicted
        5. Protein uncertain
        """
        return int(rec["proteinExistence"].split(":")[0])

    def parse_db(self, rec: dict) -> str:
        entry_type = rec["entryType"]
        if "Swiss-Prot" in entry_type:
            db = "SP"
        elif "TrEMBL" in entry_type:
            db = "TR"
        else:
            raise Exception(f"Not SP or TR, but {entry_type}")
        return db

    def parse_species(self, rec: dict) -> str:
        return rec["organism"]["scientificName"]
