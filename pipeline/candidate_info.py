import os

from Bio.PDB import PDBParser
from fetch_unbound_data import CHAINS_SEPARATOR, sub_nan
import numpy as np

from process_unbound_data import SEQUENCES, Conformation

DB_INFO_PATH = 'db_info.csv'
DB_PATH = 'data'

SEQS = 'seqs'
ANNOTATION = 'annotation'

DOT_PDB = '.pdb'
DOT_FASTA = '.fasta'

DUPLICATES_CSV = 'duplicates.csv'

CDRS = ['CDR1', 'CDR2', 'CDR3']


class CandidateInfo:
    def __init__(self, df_row, df_gaps=None):
        self.comp_name = df_row['comp_name']
        self.candidate_type = df_row['candidate_type']
        self.candidate_id = df_row['candidate_id']

        self.pdb_id_b = df_row['pdb_id_b']
        self.resolution_b = df_row['resolution_b']
        self.resolution_method_b = df_row['resolution_method_b']
        self.ab_chain_ids_b = df_row['ab_chain_ids_b'].split(CHAINS_SEPARATOR)
        self.ag_chain_ids_b = df_row['ag_chain_ids_b'].split(CHAINS_SEPARATOR)

        self.ab_pdb_id_u = df_row['ab_pdb_id_u']
        self.ab_resolution_u = df_row['ab_resolution_u']
        self.ab_resolution_method_u = df_row['ab_resolution_method_u']
        self.ab_chain_ids_u = df_row['ab_chain_ids_u'].split(CHAINS_SEPARATOR)

        self.ag_pdb_id_u = df_row['ag_pdb_id_u']
        self.ag_resolution_u = df_row['ag_resolution_u']
        self.ag_resolution_method_u = df_row['ag_resolution_method_u']
        self.ag_chain_ids_u = df_row['ag_chain_ids_u'].split(CHAINS_SEPARATOR)

        self.ab_mismatches = df_row['ab_mismatches_cnt']
        self.ag_mismatches = df_row['ag_mismatches_cnt']
        self.small_mols_msg = sub_nan(df_row['small_molecules_message'])

        self.ab_seqs_b = []
        self.ag_seqs_b = []
        self.ab_seqs_u = []
        self.ag_seqs_u = []

        self.ab_cdrs_annotation_b = []

        self.in_between_b = 0
        self.one_side_b = 0
        self.long_b = 0
        self.total_b = 0

        self.in_between_u = 0
        self.one_side_u = 0
        self.long_u = 0
        self.total_u = 0

        if df_gaps is None:
            return

        df_gaps_b, df_gaps_u = df_gaps

        selection_b = df_gaps_b['comp_name'] == self.comp_name

        selection_u = np.logical_and(df_gaps_u['comp_name'] == self.comp_name,
                                     df_gaps_u[
                                         'candidate_id'] == self.candidate_id)

        if any(selection_b):
            df_gaps_row = df_gaps_b[selection_b].iloc[0]
            self.in_between_b = int(df_gaps_row['in_between'])
            self.one_side_b = int(df_gaps_row['one_side'])
            self.long_b = int(df_gaps_row['long'])
            self.total_b = int(df_gaps_row['total'])

        if any(selection_u):
            df_gaps_row = df_gaps_u[selection_u].iloc[0]
            self.in_between_u = int(df_gaps_row['in_between'])
            self.one_side_u = int(df_gaps_row['one_side'])
            self.long_u = int(df_gaps_row['long'])
            self.total_u = int(df_gaps_row['total'])

    def to_string(self, with_candidate_id=True):
        return ','.join([self.comp_name + (('_' + self.candidate_id) if
                                           with_candidate_id else ''),
                         self.candidate_type,
                         self.pdb_id_b.upper(),
                         self.resolution_b,
                         self.resolution_method_b,
                         ':'.join(self.ab_chain_ids_b),
                         ':'.join(self.ag_chain_ids_b),
                         self.ab_pdb_id_u,
                         self.ab_resolution_u,
                         self.ab_resolution_method_u,
                         ':'.join(self.ab_chain_ids_u),
                         self.ag_pdb_id_u,
                         self.ag_resolution_u,
                         self.ag_resolution_method_u,
                         ':'.join(self.ag_chain_ids_u),
                         self.ab_mismatches,
                         self.ag_mismatches,
                         self.small_mols_msg if self.small_mols_msg
                         else 'NA',
                         str(self.in_between_b),
                         str(self.one_side_b),
                         str(self.long_b),
                         str(self.total_b),
                         str(self.in_between_u),
                         str(self.one_side_u),
                         str(self.long_u),
                         str(self.total_u)])

    def load_ab_annotation(self, db_path):
        comp_path = os.path.join(db_path, self.comp_name)

        path_to_ann_comp = os.path.join(comp_path, ANNOTATION, self.pdb_id_b +
                                        DOT_FASTA)
        path_to_fasta_ab = os.path.join(os.path.join(comp_path, ANNOTATION),
                                        self.pdb_id_b +
                                        '_ab_b' +
                                        DOT_FASTA)

        if os.path.exists(path_to_ann_comp) and \
                (not os.path.exists(path_to_fasta_ab) or
                     len(open(path_to_fasta_ab, 'r').readlines())) < 2:
            ann = read_annotation(path_to_ann_comp)

            with open(path_to_fasta_ab, 'w') as f:
                for x in self.ab_chain_ids_b:
                    for cdr in CDRS:
                        f.writelines('>{}:{}\n{}\n'.
                                     format(x, cdr, ann[(x, cdr)]))

        ab_fasta_b = read_annotation(path_to_fasta_ab)

        for x in self.ab_chain_ids_b:
            annotation = {}

            for cdr in CDRS:
                annotation[cdr] = ab_fasta_b[(x, cdr)]

            self.ab_cdrs_annotation_b.append(annotation)

    def load_sequences(self, db_path):
        comp_path = os.path.join(db_path, self.comp_name)

        ab_fasta_b = read_fasta(
            os.path.join(os.path.join(comp_path, SEQS),
                         self.pdb_id_b +
                         '_ab_b' + DOT_FASTA))
        ag_fasta_b = read_fasta(
            os.path.join(os.path.join(comp_path, SEQS),
                         self.pdb_id_b +
                         '_ag_b' + DOT_FASTA))

        self.complex_fasta_b = {**ab_fasta_b, **ag_fasta_b}

        self.ab_fasta_u = read_fasta(
            os.path.join(os.path.join(comp_path, SEQS, str(self.candidate_id)),
                         self.pdb_id_b +
                         '_ab_u' + DOT_FASTA))
        self.ag_fasta_u = read_fasta(
            os.path.join(os.path.join(comp_path, SEQS, str(self.candidate_id)),
                         self.pdb_id_b +
                         '_ag_u' + DOT_FASTA))

        for x in self.ab_chain_ids_b:
            self.ab_seqs_b.append(ab_fasta_b[x])

        for x in self.ag_chain_ids_b:
            self.ag_seqs_b.append(ag_fasta_b[x])

        for x in self.ab_chain_ids_u:
            self.ab_seqs_u.append(self.ab_fasta_u[x])

        for x in self.ag_chain_ids_u:
            self.ag_seqs_u.append(self.ag_fasta_u[x])

    # noinspection PyAttributeOutsideInit
    def to_conformation_like(self, db_path, prev_epoch):
        pdb_parser = PDBParser()

        self.comp_path = os.path.join(db_path, self.comp_name)

        epoch_path = os.path.join(self.comp_path, prev_epoch)

        self.load_sequences(db_path)

        complex_name_b = self.pdb_id_b
        complex_path_b = os.path.join(epoch_path, complex_name_b + DOT_PDB)
        self.complex_structure_b = pdb_parser.get_structure(self.comp_name,
                                                            complex_path_b)
        self.ab_chains_b = Conformation.extract_chains(
            self.complex_structure_b, self.ab_chain_ids_b)
        self.ag_chains_b = Conformation.extract_chains(
            self.complex_structure_b, self.ag_chain_ids_b)

        candidate_path = os.path.join(epoch_path, str(self.candidate_id))

        ab_name_u = self.pdb_id_b + '_ab_u'
        ab_path_u = os.path.join(candidate_path, ab_name_u + DOT_PDB)
        self.ab_structure_u = pdb_parser.get_structure('ab', ab_path_u)
        self.ab_chains_u = Conformation.extract_chains(
            self.ab_structure_u, self.ab_chain_ids_u)

        ag_name_u = self.pdb_id_b + '_ag_u'
        ag_path_u = os.path.join(candidate_path, ag_name_u + DOT_PDB)
        self.ag_structure_u = pdb_parser.get_structure('ag', ag_path_u)
        self.ag_chains_u = Conformation.extract_chains(
            self.ag_structure_u, self.ag_chain_ids_u)

        return self


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
