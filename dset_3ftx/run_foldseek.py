import argparse
import subprocess
import logging

class CommandExecutionError(Exception):
    pass

class FoldseekRunner:

    def __init__(self, query_dir, target_dir, query_db, target_db, alignment_dir, alignment_file, tmp_dir):
        self.query_dir = query_dir
        self.target_dir = target_dir
        self.query_db = query_db
        self.target_db = target_db
        self.alignment_dir = alignment_dir
        self.alignment_file = alignment_file
        self.tmp_dir = tmp_dir
        self._setup_logger()

    def __str__(self):
        return f"FoldseekRunner(query_dir={self.query_dir}, target_dir={self.target_dir}, query_db={self.query_db}, target_db={self.target_db}, alignment_dir={self.alignment_dir}, alignment_file={self.alignment_file}, tmp_dir={self.tmp_dir})"

    def _setup_logger(self):
        logging.basicConfig(level=logging.INFO)

    @classmethod
    def from_args(cls):
        parser = argparse.ArgumentParser(description='Run Foldseek commands.')
        parser.add_argument('--query_dir', required=True, help='Query directory')
        parser.add_argument('--target_dir', required=True, help='Target directory')
        parser.add_argument('--query_db', required=True, help='Query database')
        parser.add_argument('--target_db', required=True, help='Target database')
        parser.add_argument('--alignment_dir', required=True, help='Alignment directory')
        parser.add_argument('--alignment_file', required=True, help='Alignment file')
        parser.add_argument('--tmp_dir', required=True, help='Temporary directory')

        args = parser.parse_args()
        return cls(args.query_dir, args.target_dir, args.query_db, args.target_db, args.alignment_dir, args.alignment_file, args.tmp_dir)

    def run_command(self, command):
        try:
            logging.info(f"Running command: {command}")
            process = subprocess.run(command, check=True, text=True, capture_output=True, shell=True)
            logging.info(f"Command output: {process.stdout}")
        except subprocess.CalledProcessError as e:
            logging.error(f"Command '{command}' failed with error: {e.stderr}")
            raise CommandExecutionError(f"Command '{command}' failed with error: {e.stderr}")

    def run(self):
        commands = [
            f"foldseek createdb {self.query_dir} {self.query_db}",
            f"foldseek createdb {self.target_dir} {self.target_db}",
            f"foldseek search {self.query_db} {self.target_db} {self.alignment_dir} {self.tmp_dir} -a",
            f"foldseek convertalis {self.query_db} {self.target_db} {self.alignment_dir} {self.alignment_file} --format-mode 4 --format-output query,target,gapopen,pident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,cigar,qseq,tseq,qaln,taln,mismatch,qcov,tcov,lddt,qtmscore,ttmscore,alntmscore,rmsd"
        ]

        for command in commands:
            self.run_command(command)

if __name__ == "__main__":
    runner = FoldseekRunner.from_args()
