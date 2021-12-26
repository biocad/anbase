import os
import shutil

import pandas as pd
from collections import defaultdict

from candidate_info import CandidateInfo
from process_unbound_data import Conformation, ALIGNED, HETATMS_DELETED, \
    SEQUENCES
from prepper import PREPPED

ANBASE_SUMMARY_CSV = 'anbase_summary.csv'

DB_PATH = 'data'
DB_INFO_PATH = 'db_info_0.csv'
DUPLICATES_PATH = 'duplicates.csv'
GAPS_B_PATH = 'gap_stats_b.csv'
GAPS_U_PATH = 'gap_stats_u.csv'

ANBASE_SUMMARY_COLUMNS = ['comp_name', 'type',

                         'pdb_id_b', 'resolution_b', 'resolution_method_b',
                         'ab_chain_ids_b', 'ag_chain_ids_b',

                         'ab_pdb_id_u', 'ab_resolution_u',
                         'ab_resolution_method_u', 'ab_chain_ids_u',

                         'ag_pdb_id_u', 'ag_resolution_u',
                         'ag_resolution_method_u', 'ag_chain_ids_u',

                         'ab_mismatches_cnt', 'ag_mismatches_cnt',
                         'small_molecules_message', 'in_between_gaps_b',
                         'one_side_gaps_b', 'long_gaps_b', 'total_gaps_b',
                         'in_between_gaps_u',
                         'one_side_gaps_u', 'long_gaps_u', 'total_gaps_u',
                         'is_perfect']

ANBASE_SUMMARY_HEADER = ','.join(ANBASE_SUMMARY_COLUMNS)
ALTERNATIVE_CANDIDATES_COLUMNS = ['candidate_name'] + ANBASE_SUMMARY_COLUMNS[
                                                      1:-1]
ALTERNATIVE_CANDIDATES_HEADER = ','.join(ALTERNATIVE_CANDIDATES_COLUMNS)

ALTERNATIVE_CANDIDATES = 'alternative_candidates'

ANBASE_DATA_PATH = 'anbase'


def check_perfect_candidate(candidate):
    return candidate.in_between_u == 0 and candidate.one_side_u == 0 \
           and candidate.small_mols_msg is None and os.path.exists(
        os.path.join(DB_PATH, candidate.comp_name, PREPPED,
                     candidate.candidate_id))

def finalize_complex(comp_name, candidate_infos, duplicates):
    comp_candidates = candidate_infos[comp_name]
    comp_duplicates_candidates = [x for l in list(
       map(lambda x: candidate_infos[x], duplicates[comp_name])) for x in l]

    all_candidates = comp_candidates + comp_duplicates_candidates

    all_candidates.sort(key=lambda x: 0 if x.small_mols_msg is None else (
        1 if x.small_mols_msg == Conformation.MOLS_WARNING else 2))
    all_candidates.sort(key=lambda x: x.ab_mismatches + x.ag_mismatches)
    all_candidates.sort(key=lambda x: x.one_side_u)
    all_candidates.sort(key=lambda x: x.in_between_u)
    all_candidates.sort(key=lambda x: 0 if
    os.path.exists(os.path.join(DB_PATH, x.comp_name, PREPPED,
                                x.candidate_id)) else 1)

    # final_candidates = []
    # for candidate in all_candidates:
    #   is_dubplicate = any(map(
    #     lambda c: c.ab_chain_ids_u == candidate.ab_chain_ids_u and c.ag_chain_ids_u == candidate.ag_chain_ids_u, 
    #     final_candidates))
      
    #   if not is_dubplicate:
    #     final_candidates.append(candidate)


    # assign candidate idx to every candidate and
    # assign final_complex_name to put 'alternative' candidates into a single complex folder later
    i = 0
    # for candidate in all_candidates:
      # candidate.final_comp_name = final_candidates[0].comp_name
      # candidate.ab_chain_ids_b  = final_candidates[0].ab_chain_ids_b
      # candidate.ag_chain_ids_b  = final_candidates[0].ag_chain_ids_b
      
      # candidate.final_candidate_id = str(i)

      # i+=1

    is_perfect_flags = list(map(lambda x: check_perfect_candidate(x), all_candidates))

    return all_candidates, is_perfect_flags

def move_candidate_to_dir(db_path, candidate_info, dir_path):
    comp_path = os.path.join(db_path, candidate_info.comp_name)

    def move_to_dir_path(folder_name):
        path_to_folder = os.path.join(comp_path, folder_name, candidate_info.candidate_id)

        path_to_dst_unbound = os.path.join(dir_path, 
                                           folder_name.replace('prepared', 'prepared_schrod'), 
                                           candidate_info.candidate_id)

        path_to_dst_bound = os.path.join(dir_path, 
                                         folder_name.replace('prepared', 'prepared_schrod'))

        try:
            if not os.path.exists(path_to_dst_unbound):
                os.makedirs(path_to_dst_unbound)

            for file in os.listdir(path_to_folder):
                shutil.copyfile(os.path.join(path_to_folder, file),
                                os.path.join(path_to_dst_unbound, file.replace('_l_', '_ag_').replace('_r_', '_ab_')))

            if candidate_info.candidate_id == '0': # copy bound strctures only once 
                                                   # with the first pair unbound structures
              files_in_folder = filter(lambda x: candidate_info.pdb_id_b in x,
                                      os.listdir(os.path.join(comp_path, folder_name)))

              for file in files_in_folder:
                  shutil.copyfile(
                      os.path.join(comp_path, folder_name, file),
                      os.path.join(path_to_dst_bound, file))
        except Exception as e:
            print('Unsuccessful move:', path_to_folder, e, flush=True)
            return False

        return True

    move_to_dir_path(ALIGNED)
    move_to_dir_path(HETATMS_DELETED)
    prepped_moved = move_to_dir_path(PREPPED)
    move_to_dir_path(SEQUENCES)
    # move_to_dir_path(ANNOTATION)
    # move_to_dir_path(EPITOPE)

    return prepped_moved


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
    parser.add_option('--anbase-data', default=ANBASE_DATA_PATH,
                      dest='anbase_data',
                      metavar='ABASE_DATA_PATH',
                      help='Path where anbase\'s data will be stored '
                           '[default: {}]'.
                      format(DB_INFO_PATH))
    parser.add_option('--duplicates', default=DUPLICATES_PATH,
                      dest='duplicates',
                      metavar='DUPLICATES',
                      help='Path to csv with duplicates info [default: {}]'.
                      format(DUPLICATES_PATH))
    parser.add_option('--gaps-b', default=DUPLICATES_PATH,
                      dest='gaps_b',
                      metavar='GAPS_B',
                      help='Path to csv with gaps info on bound complexes'
                           ' [default: {}]'.
                      format(GAPS_B_PATH))
    parser.add_option('--gaps-u', default=DUPLICATES_PATH,
                      dest='gaps_u',
                      metavar='GAPS_B',
                      help='Path to csv with gaps info on unbound complexes '
                           '[default: {}]'.
                      format(GAPS_U_PATH))
    parser.add_option('--only-uu', default=False,
                      dest='only_uu', metavar='ONLY_UU',
                      help='Flag to process only candidates of type UU. '
                           '[default: False]')
    options, _ = parser.parse_args()

    db_df = pd.read_csv(options.db_info, dtype=str)

    gaps = None

    if os.path.exists(options.gaps_b) and os.path.exists(options.gaps_u):
        gaps_b_df = pd.read_csv(options.gaps_b, dtype=str)
        gaps_u_df = pd.read_csv(options.gaps_u, dtype=str)

    complexes = set()
    candidate_infos = defaultdict(list)

    for i in range(len(db_df)):
        candidate_info = CandidateInfo(db_df.iloc[i], gaps)

        if options.only_uu and candidate_info.candidate_type != 'U:U':
            continue

        complexes.add(candidate_info.comp_name)
        candidate_infos[candidate_info.comp_name].append(candidate_info)

    deleted = set()

    duplicates = defaultdict(list)
    dup_df = pd.read_csv(options.duplicates, dtype=str)

    for i in range(len(dup_df)):
        duplicates[dup_df.iloc[i]['comp_name']].append(dup_df.iloc[i]['duplicate_name'])

    for comp_name, duplicates_ in duplicates.items():
        if comp_name in deleted:
            continue

        for x in duplicates_:
            deleted.add(x)

            if x in complexes:
                complexes.remove(x)

    if not os.path.exists(options.anbase_data):
        os.makedirs(options.anbase_data)

    with open(ANBASE_SUMMARY_CSV, 'w') as anbase_summary_csv:
        anbase_summary_csv.write(ANBASE_SUMMARY_HEADER + '\n')
        anbase_summary_csv.flush()

        complexes_l = list(complexes)
        complexes_l.sort()

        for comp in complexes_l:
          all_candidates, is_perfect_flags = finalize_complex(comp, candidate_infos, duplicates)

          for (candidate, is_perfect) in zip(all_candidates, is_perfect_flags):

            comp_path = os.path.join(options.anbase_data, candidate.comp_name)

            successful_move = move_candidate_to_dir(options.db, candidate, comp_path)

            if successful_move:
                anbase_summary_csv.write(
                    candidate.to_string(with_candidate_id=False) +
                    ',' + str(is_perfect) + '\n')
                anbase_summary_csv.flush()
