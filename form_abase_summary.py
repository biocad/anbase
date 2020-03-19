import os
import shutil

import pandas as pd
from collections import defaultdict

from candidate_info import CandidateInfo, ANNOTATION
from collect_db_final import Conformation, ALIGNED, HETATMS_DELETED, SEQUENCES
from prepper import PREPPED, SCHROD

ABASE_SUMMARY_CSV = 'abase_summary.csv'

DB_PATH = 'data'
DB_INFO_PATH = 'db_info.csv'
DUPLICATES_PATH = 'duplicates.csv'
GAPS_PATH = 'gap_stats_u.csv'

ABASE_SUMMARY_COLUMNS = ['comp_name', 'candidate_type', 'pdb_id_b',
                         'ab_chain_ids_b', 'ag_chain_ids_b', 'ab_pdb_id_u',
                         'ab_chain_ids_u', 'ag_pdb_id_u', 'ag_chain_ids_u',
                         'ab_mismatches_cnt', 'ag_mismatches_cnt',
                         'small_molecules_message', 'in_between_gaps',
                         'one_side_gaps', 'long_gaps', 'total_gaps',
                         'is_perfect']
ABASE_SUMMARY_HEADER = ','.join(ABASE_SUMMARY_COLUMNS)

ALTERNATIVE_CANDIDATES = 'alternative_candidates'

ABASE_DATA_PATH = 'abase'


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
    all_candidates.sort(key=lambda x: x.ab_mismatches + x.ag_mismatches)
    all_candidates.sort(key=lambda x: x.one_side)
    all_candidates.sort(key=lambda x: x.in_between)

    return all_candidates[0], True, all_candidates[1:]


def move_candidate_to_dir(db_path, candidate_info, dir_path):
    comp_path = os.path.join(db_path, candidate_info.comp_name)

    def move_to_dir_path(folder_name):
        path_to_folder = os.path.join(comp_path, folder_name,
                                      candidate_info.candidate_id)

        path_to_dst = os.path.join(dir_path, folder_name)

        try:
            if not os.path.exists(path_to_dst):
                os.makedirs(path_to_dst)

            for file in os.listdir(path_to_folder):
                shutil.copyfile(os.path.join(path_to_folder, file),
                                os.path.join(path_to_dst, file))

            file_in_folder = next(
                filter(lambda x: candidate_info.pdb_id_b in x,
                       os.listdir(os.path.join(comp_path,
                                               folder_name))))

            shutil.copyfile(
                os.path.join(comp_path, folder_name, file_in_folder),
                os.path.join(path_to_dst, file_in_folder))
        except Exception as e:
            print('Unsuccessful move:', path_to_folder, e, flush=True)

    move_to_dir_path(ALIGNED)
    move_to_dir_path(HETATMS_DELETED)
    move_to_dir_path(PREPPED + '_' + SCHROD)
    move_to_dir_path(SEQUENCES)
    move_to_dir_path(ANNOTATION)


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
    parser.add_option('--abase-data', default=ABASE_DATA_PATH,
                      dest='abase_data',
                      metavar='ABASE_DATA_PATH',
                      help='Path where ABase\'s data will be stored '
                           '[default: {}]'.
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

    if not os.path.exists(options.abase_data):
        os.makedirs(options.abase_data)

    with open(ABASE_SUMMARY_CSV, 'w') as abase_summary_csv:
        abase_summary_csv.write(ABASE_SUMMARY_HEADER + '\n')
        abase_summary_csv.flush()

        for comp in complexes:
            final_candidate, is_perfect, alternative_candidates = \
                finalize_complex(comp, candidate_infos, duplicates)

            abase_summary_csv.write(
                final_candidate.to_string(with_candidate_id=False) +
                ','.join([str(is_perfect)]) + '\n')
            abase_summary_csv.flush()

            comp_path = os.path.join(options.abase_data,
                                     final_candidate.comp_name)

            move_candidate_to_dir(options.db, final_candidate, comp_path)

            with open(os.path.join(comp_path, ALTERNATIVE_CANDIDATES + '.csv'),
                      'w') as alternative_candidates_csv:
                alternative_candidates_csv.write(
                    ','.join(ABASE_SUMMARY_COLUMNS[:-2]) + '\n')

                if alternative_candidates is not None:
                    for alternative_candidate in alternative_candidates:
                        alternative_candidates_csv.write(
                            alternative_candidate.to_string() + '\n')

            if alternative_candidates is None:
                continue

            alternative_candidates_path = os.path.join(comp_path,
                                                       ALTERNATIVE_CANDIDATES)

            if not os.path.exists(alternative_candidates_path):
                os.mkdir(alternative_candidates_path)

            for alternative_candidate in alternative_candidates:
                path_to_alternative_candidate = os.path.join(
                    alternative_candidates_path,
                    alternative_candidate.comp_name)

                if not os.path.exists(path_to_alternative_candidate):
                    os.mkdir(path_to_alternative_candidate)

                move_candidate_to_dir(options.db, alternative_candidate,
                                      path_to_alternative_candidate)
