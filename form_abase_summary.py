import pandas as pd
from collections import defaultdict
import math
import numpy as np

from collect_db_final import Conformation

ABASE_SUMMARY_CSV = 'abase_summary.csv'

DB_PATH = 'data'
DB_INFO_PATH = 'db_info.csv'
DUPLICATES_PATH = 'duplicates.csv'
GAPS_PATH = 'gap_stats_u.csv'


def sub_nan(val):
    if isinstance(val, float) and math.isnan(val):
        return None
    return val


class CandidateInfo:
    def __init__(self, df_row, df_gaps):
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

        self.small_mols_msg = sub_nan(df_row['small_molecules_message'])

        selection = np.logical_and(df_gaps['comp_name'] == self.comp_name,
                                   df_gaps[
                                       'candidate_id'] == self.candidate_id)

        self.in_between = 0
        self.one_side = 0
        self.long = 0
        self.total = 0

        if any(selection):
            df_gaps_row = df_gaps[selection].iloc[0]
            self.in_between = int(df_gaps_row['in_between'])
            self.one_side = int(df_gaps_row['one_side'])
            self.long = int(df_gaps_row['long'])
            self.total = int(df_gaps_row['total'])

    def to_string(self, with_candidate_id=True):
        addition = [self.candidate_id] if with_candidate_id else []
        return ','.join([self.comp_name.replace(':', '+'), self.candidate_type] + addition +
                        [self.pdb_id_b.upper(),
                         ':'.join(self.ab_chain_ids_b),
                         ':'.join(self.ag_chain_ids_b),
                         self.ab_pdb_id_u,
                         ':'.join(self.ab_chain_ids_u),
                         self.ag_pdb_id_u,
                         ':'.join(self.ag_chain_ids_u),
                         self.small_mols_msg if self.small_mols_msg
                         else 'NA',
                         str(self.in_between),
                         str(self.one_side),
                         str(self.long),
                         str(self.total)])


def check_perfect_candidate(candidate):
    return candidate.in_between == 0 and candidate.one_side == 0 \
           and candidate.small_mols_msg is None


def finalize_complex(comp_name, candidate_infos, duplicates):
    comp_candidates = candidate_infos[comp_name]
    comp_duplicates_candidates = [x for l in list(
        map(lambda x: candidate_infos[x], duplicates[comp_name])) for x in l]

    all_candidates = comp_candidates + comp_duplicates_candidates

    ideal_candidates = list(
        filter(lambda x: check_perfect_candidate(x), all_candidates))

    if len(ideal_candidates) > 0:
        return ideal_candidates[0], True, None

    all_candidates.sort(key=lambda x: 0 if x.small_mols_msg is None else (
        1 if x.small_mols_msg == Conformation.MOLS_WARNING else 2))
    all_candidates.sort(key=lambda x: x.one_side)
    all_candidates.sort(key=lambda x: x.in_between)

    return all_candidates[0], True, all_candidates[1:]


if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('--db', default=DB_PATH, dest='db', metavar='DB',
                      help='Path to dev database [default: {}]'.format(
                          DB_PATH))
    parser.add_option('--db-info', default=DB_INFO_PATH, dest='db_info',
                      metavar='DB_INFO_PATH',
                      help='Path to dev database info csv file [default: {}]'.
                      format(DB_INFO_PATH))
    parser.add_option('--duplicates', default=DUPLICATES_PATH,
                      dest='duplicates',
                      metavar='DUPLICATES',
                      help='Path to csv with duplicates info [default: {}]'.
                      format(DUPLICATES_PATH))
    parser.add_option('--gaps', default=DUPLICATES_PATH,
                      dest='gaps',
                      metavar='GAPS',
                      help='Path to csv with gaps info [default: {}]'.
                      format(GAPS_PATH))
    parser.add_option('--only-uu', default=False,
                      dest='only_uu', metavar='ONLY_UU',
                      help='Flag to process only candidates of type UU. '
                           '[default: False]')
    options, _ = parser.parse_args()

    db_df = pd.read_csv(options.db_info, dtype=str)
    gaps_df = pd.read_csv(options.gaps, dtype=str)

    complexes = set()
    candidate_infos = defaultdict(list)

    for i in range(len(db_df)):
        candidate_info = CandidateInfo(db_df.iloc[i], gaps_df)

        if options.only_uu and candidate_info.candidate_type != 'U:U':
            continue

        complexes.add(candidate_info.comp_name)
        candidate_infos[candidate_info.comp_name].append(candidate_info)

    deleted = set()

    duplicates = defaultdict(list)
    dup_df = pd.read_csv(options.duplicates, dtype=str)

    for i in range(len(dup_df)):
        duplicates[dup_df.iloc[i]['comp_name']].append(
            dup_df.iloc[i]['duplicate_name'])

    for comp_name, duplicates_ in duplicates.items():
        if comp_name in deleted:
            continue

        for x in duplicates_:
            deleted.add(x)

            if x in complexes:
                complexes.remove(x)

    with open(ABASE_SUMMARY_CSV, 'w') as abase_summary_csv:
        abase_summary_csv.write(
            'comp_name,candidate_type,pdb_id_b,'
            'ab_chain_ids_b,ag_chain_ids_b,ab_pdb_id_u,ab_chain_ids_u,ag_'
            'pdb_id_u,ag_chain_ids_u,small_molecules_message,'
            'in_between_gaps,one_side_gaps,long_gaps,total_gaps\n')
        abase_summary_csv.flush()

        for comp in complexes:
            final_candidate, is_perfect, second_choices = \
                finalize_complex(comp, candidate_infos, duplicates)

            abase_summary_csv.write(
                final_candidate.to_string(with_candidate_id=False) + '\n')
            abase_summary_csv.flush()
