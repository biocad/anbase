import os

from collect_db import CHAINS_SEPARATOR, sub_nan
import numpy as np

DB_INFO_PATH = 'db_info.csv'
DB_PATH = 'data'

SEQS = 'seqs'
ANNOTATION = 'annotation'

DOT_FASTA = '.fasta'

DUPLICATES_CSV = 'duplicates.csv'

CDRS = ['CDR1', 'CDR2', 'CDR3']


class CandidateInfo:
    def __init__(self, df_row, df_gaps=None):
        self.comp_name = df_row['comp_name']
        self.candidate_type = df_row['candidate_type']
        self.candidate_id = df_row['candidate_id']

        self.pdb_id_b = df_row['pdb_id_b']
        self.ab_chain_ids_b = df_row['ab_chain_ids_b'].split(CHAINS_SEPARATOR)
        self.ag_chain_ids_b = df_row['ag_chain_ids_b'].split(CHAINS_SEPARATOR)

        self.ab_pdb_id_u = df_row['ab_pdb_id_u']
        self.ab_chain_ids_u = df_row['ab_chain_ids_u'].split(CHAINS_SEPARATOR)
        self.ag_pdb_id_u = df_row['ag_pdb_id_u']
        self.ag_chain_ids_u = df_row['ag_chain_ids_u'].split(CHAINS_SEPARATOR)

        self.ab_mismatches = df_row['ab_mismatches_cnt']
        self.ag_mismatches = df_row['ag_mismatches_cnt']
        self.small_mols_msg = sub_nan(df_row['small_molecules_message'])

        self.ab_seqs = []
        self.ag_seqs = []

        self.ab_cdrs_annotation_b = []

        self.in_between = 0
        self.one_side = 0
        self.long = 0
        self.total = 0

        if df_gaps is None:
            return

        selection = np.logical_and(df_gaps['comp_name'] == self.comp_name,
                                   df_gaps[
                                       'candidate_id'] == self.candidate_id)

        if any(selection):
            df_gaps_row = df_gaps[selection].iloc[0]
            self.in_between = int(df_gaps_row['in_between'])
            self.one_side = int(df_gaps_row['one_side'])
            self.long = int(df_gaps_row['long'])
            self.total = int(df_gaps_row['total'])

    def to_string(self, with_candidate_id=True):
        addition = [self.candidate_id] if with_candidate_id else []
        return ','.join([self.comp_name.replace(':', '+'),
                         self.candidate_type] + addition +
                        [self.pdb_id_b.upper(),
                         ':'.join(self.ab_chain_ids_b),
                         ':'.join(self.ag_chain_ids_b),
                         self.ab_pdb_id_u,
                         ':'.join(self.ab_chain_ids_u),
                         self.ag_pdb_id_u,
                         ':'.join(self.ag_chain_ids_u),
                         self.ab_mismatches,
                         self.ag_mismatches,
                         self.small_mols_msg if self.small_mols_msg
                         else 'NA',
                         str(self.in_between),
                         str(self.one_side),
                         str(self.long),
                         str(self.total)])

    def load_ab_annotation(self, db_path):
        comp_path = os.path.join(db_path, self.comp_name)
        ab_fasta_b = read_annotation(
            os.path.join(os.path.join(comp_path, ANNOTATION), self.pdb_id_b +
                         DOT_FASTA))

        for x in self.ab_chain_ids_b:
            annotation = {}

            for cdr in CDRS:
                annotation[cdr] = ab_fasta_b[(x, cdr)]

            self.ab_cdrs_annotation_b.append(annotation)

    def load_sequences(self, db_path):
        comp_path = os.path.join(db_path, self.comp_name)

        complex_fasta_b = read_fasta(
            os.path.join(os.path.join(comp_path, SEQS), self.pdb_id_b +
                         DOT_FASTA))

        for x in self.ab_chain_ids_b:
            self.ab_seqs.append(complex_fasta_b[x])

        for x in self.ag_chain_ids_b:
            self.ag_seqs.append(complex_fasta_b[x])


def read_fasta(path):
    res = {}

    with open(path, 'r') as f:
        lines = f.readlines()

        i = 0
        while i < len(lines):
            res[lines[i].split(':')[1].strip()] = lines[i + 1].strip()
            i += 2

    return res


def read_annotation(path):
    res = {}

    with open(path, 'r') as f:
        lines = f.readlines()

        i = 0
        while i < len(lines):
            [chain_id, region] = lines[i].strip()[1:].split(':')
            res[(chain_id, region)] = lines[i + 1].strip()
            i += 2

    return res