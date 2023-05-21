#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on: Tue 09 May 2023 18:29:23
Description: Download best PDB structure from a having UniProt JSON files
Usage:       python dset_3ftx/pdb_get_structures.py data/uniprot_entries/ data/pdb --verbose

@author: tsenoner
"""
import argparse
import json
import re
from pathlib import Path

import requests


class PdbEntry:
    def __init__(self, entry_id, method=None, resolution=None, position=None):
        self.entry_id = entry_id
        self.method = method
        self.resolution = resolution
        self.position = position
        self.sequence = None

    def set_sequence(self, sequence):
        self.sequence = sequence[int(self.position[1]) : int(self.position[2])]

    def __getitem__(self, item):
        return getattr(self, item)

    def __repr__(self):
        return (
            f"<PdbEntry id={self.entry_id}, method={self.method},"
            f" resolution={self.resolution}, position={self.position},"
            f" sequence={self.sequence}>"
        )


class PdbData:
    def __init__(self, uniprot_dir, verbose=False):
        self.uniprot_dir = uniprot_dir
        self.verbose = verbose
        self.data = self.process_uniprot_data()

    @staticmethod
    def read_json_file(json_file: str):
        """Read a JSON file and return the data."""
        with open(json_file, "r") as json_handle:
            raw_data = json.load(json_handle)
        return raw_data

    @staticmethod
    def extract_pdb_entries(cross_references):
        """Extract PDB entries from UniProt cross references."""

        def parse_properties(properties):
            method, resolution, position = None, None, None
            for prop in properties:
                if prop["key"] == "Method":
                    method = prop["value"]
                if prop["key"] == "Resolution" and prop["value"] != "-":
                    resolution = float(prop["value"].split(" ")[0])
                if prop["key"] == "Chains":
                    position = re.match(
                        r"([\w/]+)=(\d+)-(\d+)", prop["value"]
                    ).groups()
            return method, resolution, position

        return [
            PdbEntry(
                cross_ref["id"], *parse_properties(cross_ref["properties"])
            )
            for cross_ref in cross_references
            if cross_ref["database"] == "PDB"
        ]

    @staticmethod
    def choose_best_pdb_entry(pdb_entries, sequence):
        """Choose the best PDB entry based on resolution."""
        if not pdb_entries:
            return None

        for entry in pdb_entries:
            entry.set_sequence(sequence)

        # filter for longest sequence
        max_len = len(max(pdb_entries, key=lambda x: len(x.sequence)).sequence)
        pdb_entries = [entry for entry in pdb_entries if len(entry.sequence) == max_len]

        return min(
            pdb_entries, key=lambda x: (x.resolution is None, x.resolution)
        )

    def process_uniprot_data(self):
        data = {}
        for json_file in Path(self.uniprot_dir).glob("*.json"):
            raw_data = self.read_json_file(json_file=json_file)

            if "uniProtKBCrossReferences" not in raw_data:
                continue

            uniprot_id = raw_data["primaryAccession"]
            sequence = raw_data["sequence"]["value"]
            cross_references = raw_data["uniProtKBCrossReferences"]

            pdb_entries = self.extract_pdb_entries(cross_references)
            best_entry = self.choose_best_pdb_entry(pdb_entries, sequence)

            if best_entry is not None:
                data[uniprot_id] = best_entry
        return data

    def download_pdb_file(self, pdb_id, output_dir="pdb_files"):
        """Download a PDB file from the RCSB PDB server and save it in the output directory."""
        output_path = Path(output_dir) / f"{pdb_id}.pdb"

        if output_path.is_file():
            if self.verbose:
                print(
                    f"PDB file for ID: {pdb_id} already exists. Skipping"
                    " download."
                )
            return

        pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        response = requests.get(pdb_url)
        if response.status_code == 200:
            # Ensure the output directory exists
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            # Save the PDB file
            with open(output_path, "w") as pdb_file:
                pdb_file.write(response.text)
            if self.verbose:
                print(f"Downloaded PDB file for ID: {pdb_id}")
        else:
            if self.verbose:
                print(f"Failed to download PDB file for ID: {pdb_id}")


def download_pdb_from_uniprot(uniprot_directory, pdb_output_dir, verbose=False):
    pdb_data = PdbData(uniprot_directory, verbose=verbose)
    data = pdb_data.data

    # Download PDB files for the extracted PDB ids
    for pdb_entry in data.values():
        pdb_data.download_pdb_file(
            pdb_entry.entry_id, output_dir=pdb_output_dir
        )

    return pdb_data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Download best PDB structure from UniProt JSON files"
    )
    parser.add_argument(
        "uniprot_directory",
        type=str,
        help="Path to the directory containing the UniProt JSON entries",
    )
    parser.add_argument(
        "pdb_output_dir",
        type=str,
        help=(
            "Path to the output directory where the PDB files will be"
            " downloaded"
        ),
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print download progress messages (optional)",
    )

    args = parser.parse_args()

    download_pdb_from_uniprot(
        args.uniprot_directory, args.pdb_output_dir, args.verbose
    )
