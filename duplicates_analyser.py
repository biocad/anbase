import pandas as pd
import os

from Bio import pairwise2

from collect_db_final import comp_name_to_dir_name

DB_INFO_PATH = 'db_info.csv'
DB_PATH = 'data'

SEQS = 'seqs'
ANNOTATION = 'annotation'

DOT_FASTA = '.fasta'

DUPLICATES_CSV = 'duplicates.csv'

CDRS = ['CDR1', 'CDR2', 'CDR3']


class CandidateInfo:
    def __init__(self, df_row):
        self.comp_name = df_row['comp_name']
        self.candidate_type = df_row['candidate_type']
        self.candidate_id = df_row['candidate_id']

        self.pdb_id_b = df_row['pdb_id_b']
        self.ab_chain_ids_b = df_row['ab_chain_ids_b'].split(':')
        self.ag_chain_ids_b = df_row['ag_chain_ids_b'].split(':')

        self.ab_pdb_id_u = df_row['ab_pdb_id_u']
        self.ab_chain_ids_u = df_row['ab_chain_ids_u'].split(':')
        self.ag_pdb_id_u = df_row['ag_pdb_id_u']
        self.ag_chain_ids_u = df_row['ag_chain_ids_u'].split(':')

        self.ab_seqs = []
        self.ag_seqs = []

        self.ab_cdrs_annotation_b = []

    def load_ab_annotation(self, db_path):
        comp_path = os.path.join(db_path,
                                 comp_name_to_dir_name(self.comp_name))
        ab_fasta_b = read_annotation(
            os.path.join(os.path.join(comp_path, ANNOTATION), self.pdb_id_b +
                         DOT_FASTA))

        for x in self.ab_chain_ids_b:
            annotation = {}

            for cdr in CDRS:
                annotation[cdr] = ab_fasta_b[(x, cdr)]

            self.ab_cdrs_annotation_b.append(annotation)

    def load_sequences(self, db_path):
        comp_path = os.path.join(db_path,
                                 comp_name_to_dir_name(self.comp_name))

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


def calc_matches_mismatches(seq1, seq2):
    alignment_list = pairwise2.align.localxs(seq1, seq2, -1, -1,
                                             penalize_end_gaps=False,
                                             one_alignment_only=True)

    if not alignment_list:
        return False

    alignment = alignment_list[0]

    matches_count = 0

    query_alignment = alignment[0]
    target_alignment = alignment[1]

    for i in range(len(query_alignment)):
        if query_alignment[i] == target_alignment[i]:
            matches_count += 1

    return matches_count, len(query_alignment) - matches_count


def similarity_of_two_seqs(seq1, seq2):
    matches, mismatches = calc_matches_mismatches(seq1, seq2)

    score = float(matches) / float(matches + mismatches)

    return score >= 0.9


def similarity_of_abs(comp1, comp2):
    for i in range(len(comp1.ab_chain_ids_b)):
        for cdr in CDRS:
            seq1 = comp1.ab_cdrs_annotation_b[i][cdr]
            seq2 = comp2.ab_cdrs_annotation_b[i][cdr]

            _, mismatches = calc_matches_mismatches(seq1, seq2)

            if mismatches >= 2:
                print('Not equal:', seq1, 'and', seq2, ', mismatches:', mismatches,
                      flush=True)
                return False
    return True


def similarity_of_two_complexes(comp1, comp2):
    if len(comp1.ab_seqs) != len(comp2.ab_seqs) or \
            len(comp1.ag_seqs) != len(comp2.ag_seqs):
        return False

    ab_chains_similar = similarity_of_abs(comp1, comp2)
    ag_chains_similar = all(map(lambda p: similarity_of_two_seqs(p[0], p[1]),
                                zip(comp1.ag_seqs, comp2.ag_seqs)))

    return ab_chains_similar and ag_chains_similar


if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('--db', default=DB_PATH, dest='db', metavar='DB',
                      help='Path to database [default: {}]'.format(DB_PATH))
    parser.add_option('--db-info', default=DB_INFO_PATH, dest='db_info',
                      metavar='DB_INFO_PATH',
                      help='Path to database info csv file [default: {}]'.
                      format(DB_INFO_PATH))
    parser.add_option('--only-uu', default=False,
                      dest='only_uu', metavar='ONLY_UU',
                      help='Flag to process only candidates of type UU. '
                           '[default: False]')
    options, _ = parser.parse_args()

    df = pd.read_csv(options.db_info, dtype=str)

    complexes = set()
    complexes_with_chains = []

    for i in range(len(df)):
        candidate_info = CandidateInfo(df.iloc[i])

        if options.only_uu and candidate_info.candidate_type != 'U:U':
            continue

        if candidate_info.comp_name in complexes:
            continue

        complexes.add(candidate_info.comp_name)

        candidate_info.load_sequences(options.db)
        candidate_info.load_ab_annotation(options.db)
        complexes_with_chains.append(candidate_info)

    with open(DUPLICATES_CSV, 'w') as duplicates_csv:
        duplicates_csv.write('comp_name,duplicate_name\n')
        duplicates_csv.flush()

        for comp in complexes_with_chains:
            similar_comps = []
            for other_comp in complexes_with_chains:
                if other_comp.comp_name == comp.comp_name:
                    continue

                if similarity_of_two_complexes(comp, other_comp):
                    similar_comps.append(other_comp.comp_name)

            for similar_comp in similar_comps:
                duplicates_csv.write(
                    '{},{}\n'.format(comp.comp_name, similar_comp))
                duplicates_csv.flush()
