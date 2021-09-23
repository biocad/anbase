from multiprocessing.pool import Pool

import pandas as pd

from candidate_info import CandidateInfo
from fetch_unbound_data import is_subsequence_of

DB_INFO_PATH = 'db_info_0.csv'
DB_PATH = 'data'

SEQS = 'seqs'
ANNOTATION = 'annotation'

DOT_FASTA = '.fasta'

DUPLICATES_CSV = 'duplicates.csv'

CDRS = ['CDR1', 'CDR2', 'CDR3']

NUMBER_OF_PROCESSES = 50


def similarity_of_two_seqs(seq1, seq2):
    return is_subsequence_of(seq1, seq2)

def similarity_of_two_complexes(comp1, comp2):
    if len(comp1.ab_seqs_b) != len(comp2.ab_seqs_b) or \
            len(comp1.ag_seqs_b) != len(comp2.ag_seqs_b):
        return False

    # TODO: check if it works
    if len(comp1.ab_seqs_b) == 2 and len(comp1.ab_seqs_b) == 2: # compare combinations of seqeunces
      comp1_fst_seq_similarity = similarity_of_two_seqs(comp1.ab_seqs_b[0], comp2.ab_seqs_b[0]) or \
                                 similarity_of_two_seqs(comp1.ab_seqs_b[0], comp2.ab_seqs_b[1])
      comp1_snd_seq_similarity = similarity_of_two_seqs(comp1.ab_seqs_b[1], comp2.ab_seqs_b[1]) or \
                                 similarity_of_two_seqs(comp1.ab_seqs_b[1], comp2.ab_seqs_b[0])
      ab_chains_similar = comp1_fst_seq_similarity and comp1_snd_seq_similarity
    else: # compare just two sequnces for single chain antibodies
      ab_chains_similar = similarity_of_two_seqs(comp1.ab_seqs_b[0], comp2.ab_seqs_b[0])

    ag_chains_similar = all(map(lambda p: similarity_of_two_seqs(p[0], p[1]),
                                zip(comp1.ag_seqs_b, comp2.ag_seqs_b)))

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

        # if candidate_info.comp_name not in ['3u2s_H:L|G', '3u4e_H:L|G']:
        #     continue

        complexes.add(candidate_info.comp_name)

        candidate_info.load_sequences(options.db)
        # candidate_info.load_ab_annotation(options.db)
        complexes_with_chains.append(candidate_info)

    with open(DUPLICATES_CSV, 'w') as duplicates_csv:
        duplicates_csv.write('comp_name,duplicate_name\n')
        duplicates_csv.flush()

        for comp in complexes_with_chains:
            print('Searching for duplicates of:', comp.comp_name, flush=True)

            similar_comps = []

            with Pool(NUMBER_OF_PROCESSES) as pool:
                res = pool.starmap(similarity_of_two_complexes,
                                      list((comp, x) for x in
                                           complexes_with_chains))

                for other_comp, are_similar in zip(complexes_with_chains, res):
                    if other_comp.comp_name == comp.comp_name:
                        continue

                    if are_similar:
                        similar_comps.append(other_comp.comp_name)

            for similar_comp in similar_comps:
                duplicates_csv.write(
                    '{},{}\n'.format(comp.comp_name, similar_comp))
                duplicates_csv.flush()
