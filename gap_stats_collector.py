import os
from collections import defaultdict

from Bio import pairwise2
from Bio.PDB import PDBParser

from candidate_info import CandidateInfo
from collect_db import comp_name_to_pdb_and_chains
from collect_db_final import SEQUENCES, HETATMS_DELETED, comp_name_to_dir_name, \
    INTERFACE_CUTOFF, extract_seq
import pandas as pd
import numpy as np

DB_PATH = 'data'
DOT_PDB = '.pdb'
DOT_FASTA = '.fasta'

DB_INFO_PATH = 'db_info.csv'

MAGIC_INTERFACE_EXTENSION_CONSTANT = 5

MAX_NOT_SUSPICIOUS_LENGTH_OF_GAP = 15

pdb_parser = PDBParser()

GAP_STATS_B_CSV = 'gap_stats_b.csv'
GAP_STATS_U_CSV = 'gap_stats_u.csv'


def get_chains(structure, chain_ids):
    return list(map(lambda x: get_chain_with_id(structure, x), chain_ids))


def union_gap_stats(gap_stats1, gap_stats2):
    return gap_stats1[0] + gap_stats2[0], \
           gap_stats1[1] + gap_stats2[1], \
           gap_stats1[2] + gap_stats2[2], \
           gap_stats1[3] + gap_stats2[3],


def get_gap_stats(complex_structure_b, ab_chain_ids_b, ag_chain_ids_b,
                  ab_structure_u,
                  ab_chain_ids_u, ag_structure_u, ag_chain_ids_u,
                  complex_fasta_b, ab_fasta_u, ag_fasta_u):
    ab_chains_b = get_chains(complex_structure_b, ab_chain_ids_b)
    ag_chains_b = get_chains(complex_structure_b, ag_chain_ids_b)
    ab_chains_u = get_chains(ab_structure_u, ab_chain_ids_u)
    ag_chains_u = get_chains(ag_structure_u, ag_chain_ids_u)

    ab_interface_residues_inds_b, ag_interface_residues_inds_b = \
        interface_residue_ids(ab_chains_b, ag_chains_b)

    ab_interface_residues_inds_u, ag_interface_residues_inds_u = \
        interface_residue_ids(list(ab_structure_u.get_chains()),
                              list(ag_structure_u.get_chains()))

    gap_stats_ab_b = gap_stats_for_chains(complex_fasta_b,
                                          ab_interface_residues_inds_b,
                                          ab_chains_b)
    gap_stats_ag_b = gap_stats_for_chains(complex_fasta_b,
                                          ag_interface_residues_inds_b,
                                          ag_chains_b)
    gap_stats_ab_u = gap_stats_for_chains(ab_fasta_u,
                                          ab_interface_residues_inds_u,
                                          ab_chains_u)
    gap_stats_ag_u = gap_stats_for_chains(ag_fasta_u,
                                          ag_interface_residues_inds_u,
                                          ag_chains_u)

    return union_gap_stats(gap_stats_ab_b, gap_stats_ag_b), union_gap_stats(
        gap_stats_ab_u, gap_stats_ag_u)


def get_gap_stats_for_chain(seq, chain, interface_residues):
    gapped_seq = extract_seq(chain)

    alignment = \
        pairwise2.align.localxs(gapped_seq, seq, -1, 0,
                                one_alignment_only=True)[0]

    # UNCOMMENT TO DEBUG
    #
    # c = -1
    # s = ''
    # for x in alignment[0]:
    #     if x != '-':
    #         c += 1
    #
    #     if x != '-' and c in interface_residues:
    #         s += '+'
    #     else:
    #         s += ' '
    #
    # print(s)
    # print(alignment[0])
    # print(alignment[1])

    gaps = []

    cur_left_bound = None
    ind = -1

    in_between = 0
    one_side = 0

    left_bound_true = 0

    for i in range(len(alignment[0])):
        symbol_gapped = alignment[0][i]

        if symbol_gapped != '-':
            ind += 1

        if symbol_gapped == '-' and cur_left_bound is None:
            cur_left_bound = ind
            left_bound_true = i
        elif symbol_gapped != '-' and cur_left_bound is not None:
            gaps.append((cur_left_bound, ind, i - left_bound_true + 1))
            cur_left_bound = None

    if cur_left_bound is not None:
        gaps.append((cur_left_bound, len(alignment[0]),
                     len(alignment[0]) - left_bound_true))

    for gap in gaps:
        if gap[0] in interface_residues and gap[1] in interface_residues:
            in_between += 1
        elif gap[0] in interface_residues or gap[1] in interface_residues:
            one_side += 1

    long_gaps = 0

    for gap in gaps:
        if gap[2] > MAX_NOT_SUSPICIOUS_LENGTH_OF_GAP:
            long_gaps += 1

    return in_between, one_side, long_gaps, len(gaps)


def get_chain_with_id(structure, chain_id):
    for chain in structure.get_chains():
        if chain.get_id() == chain_id:
            return chain


def gap_stats_for_chains(fasta, interface_residues_for_chains, chains):
    in_between_cnt = 0
    one_side_cnt = 0
    long_gaps_cnt = 0
    total_cnt = 0

    for i in range(len(chains)):
        in_between, one_side, long_gaps, total = get_gap_stats_for_chain(
            fasta[chains[i].id],
            chains[i],
            interface_residues_for_chains[chains[i].id])
        in_between_cnt += in_between
        one_side_cnt += one_side
        long_gaps_cnt += long_gaps
        total_cnt += total

    return in_between_cnt, \
           one_side_cnt, \
           long_gaps_cnt, \
           total_cnt


def interface_residue_ids(ab_chains, ag_chains):
    ab_chain_to_interface_residues = defaultdict(set)
    ag_chain_to_interface_residues = defaultdict(set)

    for ab_chain in ab_chains:
        for ag_chain in ag_chains:
            ab_ind = -1
            for ab_residue in ab_chain:
                ab_ind += 1
                ag_ind = -1
                for ag_residue in ag_chain:
                    ag_ind += 1
                    for ab_at in ab_residue:
                        if ab_at.get_id() != 'CA':
                            continue

                        can_break = False

                        for ag_at in ag_residue:
                            if ag_at.get_id() != 'CA':
                                continue

                            if np.linalg.norm(
                                    ab_at.coord - ag_at.coord) < \
                                    INTERFACE_CUTOFF + \
                                    MAGIC_INTERFACE_EXTENSION_CONSTANT:
                                ab_chain_to_interface_residues[
                                    ab_chain].add(ab_ind)
                                ag_chain_to_interface_residues[
                                    ag_chain].add(ag_ind)
                                can_break = True
                                break

                        if can_break:
                            break
    ab_res = {}
    ag_res = {}

    for ab_chain in ab_chains:
        ab_res[ab_chain.id] = ab_chain_to_interface_residues[ab_chain]

    for ag_chain in ag_chains:
        ag_res[ag_chain.id] = ag_chain_to_interface_residues[ag_chain]

    return ab_res, ag_res


def read_fasta(path):
    res = {}

    with open(path, 'r') as f:
        lines = f.readlines()

        i = 0
        while i < len(lines):
            res[lines[i].split(':')[1].strip()] = lines[i + 1].strip()
            i += 2

    return res


def process_candidate(candidate, db_path, prev_epoch, gap_stats_b_csv,
                      processed_comps,
                      gap_stats_u_csv):
    pdb_id_b, _, _ = comp_name_to_pdb_and_chains(candidate.comp_name)

    comp_path = os.path.join(db_path,
                             comp_name_to_dir_name(candidate.comp_name))

    epoch_path = os.path.join(comp_path, prev_epoch)
    seqs_path = os.path.join(comp_path, SEQUENCES)

    complex_name_b = pdb_id_b
    complex_path_b = os.path.join(epoch_path, complex_name_b + DOT_PDB)
    complex_structure_b = pdb_parser.get_structure(candidate.comp_name,
                                                   complex_path_b)
    complex_fasta_b = read_fasta(
        os.path.join(seqs_path, complex_name_b + DOT_FASTA))

    candidate_path = os.path.join(epoch_path, str(candidate.candidate_id))

    [ab_suf, ag_suf] = list(map(lambda x: '_u' if x == 'U' else '_b',
                                candidate.candidate_type.split(':')))

    ab_name_u = pdb_id_b + '_r' + ab_suf
    ab_path_u = os.path.join(candidate_path, ab_name_u + DOT_PDB)
    ab_structure_u = pdb_parser.get_structure('r', ab_path_u)
    ab_fasta_u = read_fasta(
        os.path.join(seqs_path, str(candidate.candidate_id),
                     ab_name_u + DOT_FASTA))

    ag_name_u = pdb_id_b + '_l' + ag_suf
    ag_path_u = os.path.join(candidate_path, ag_name_u + DOT_PDB)
    ag_structure_u = pdb_parser.get_structure('l', ag_path_u)
    ag_fasta_u = read_fasta(
        os.path.join(seqs_path, str(candidate.candidate_id),
                     ag_name_u + DOT_FASTA))

    gap_stats_b, gap_stats_u = get_gap_stats(complex_structure_b,
                                             candidate.ab_chain_ids_b,
                                             candidate.ag_chain_ids_b,
                                             ab_structure_u,
                                             candidate.ab_chain_ids_u,
                                             ag_structure_u,
                                             candidate.ag_chain_ids_u,
                                             complex_fasta_b, ab_fasta_u,
                                             ag_fasta_u)

    if candidate.comp_name not in processed_comps:
        processed_comps.add(candidate.comp_name)
        gap_stats_b_csv.write(
            '{},{},{},{},{}\n'.format(candidate.comp_name,
                                      gap_stats_b[0],
                                      gap_stats_b[1], gap_stats_b[2],
                                      gap_stats_b[3]))
        gap_stats_b_csv.flush()

    gap_stats_u_csv.write(
        '{},{},{},{},{},{}\n'.format(candidate.comp_name,
                                     candidate.candidate_id, gap_stats_u[0],
                                     gap_stats_u[1], gap_stats_u[2],
                                     gap_stats_u[3]))
    gap_stats_u_csv.flush()


if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('--db', default=DB_PATH, dest='db', metavar='DB',
                      help='Path to database [default: {}]'.format(DB_PATH))
    parser.add_option('--db-info', default=DB_INFO_PATH, dest='db_info',
                      metavar='DB_INFO_PATH',
                      help='Path to database info csv file [default: {}]'.
                      format(DB_INFO_PATH))
    parser.add_option('--prev-epoch', default=HETATMS_DELETED,
                      dest='prev_epoch', metavar='PREV_EPOCH',
                      help='Name of the epoch structures from which will be '
                           'checked for gaps. [default: {}]'.format(
                          HETATMS_DELETED))
    parser.add_option('--only-uu', default=False,
                      dest='only_uu', metavar='ONLY_UU',
                      help='Flag to process only candidates of type UU. '
                           '[default: False]')
    options, _ = parser.parse_args()

    processed_comps = set([])
    processed_candidates = frozenset([])

    gap_stats_b_csv_path = options.prev_epoch + '_' + GAP_STATS_B_CSV
    gap_stats_u_csv_path = options.prev_epoch + '_' + GAP_STATS_U_CSV

    if os.path.exists(gap_stats_b_csv_path):
        with open(gap_stats_b_csv_path, 'r') as f:
            lines = f.readlines()[1:]
            processed_comps = set(map(lambda x: x.split(',')[0], lines))
    else:
        header_b = 'comp_name,in_between,one_side,long,total\n'

        with open(gap_stats_b_csv_path, 'w') as f:
            f.write(header_b)
            f.flush()

    if os.path.exists(gap_stats_u_csv_path):
        with open(gap_stats_u_csv_path, 'r') as f:
            lines = f.readlines()[1:]
            processed_candidates = frozenset(
                map(lambda x: '_'.join(x.split(',')[:2]), lines))
    else:
        header_u = 'comp_name,candidate_id,in_between,one_side,long,total\n'

        with open(gap_stats_u_csv_path, 'w') as f:
            f.write(header_u)
            f.flush()

    with open(gap_stats_b_csv_path, 'a') as gap_stats_csv_b, \
            open(gap_stats_u_csv_path, 'a') as gap_stats_csv_u:

        df = pd.read_csv(options.db_info, dtype=str)

        for i in range(len(df)):
            candidate_info = CandidateInfo(df.iloc[i])

            if options.only_uu and candidate_info.candidate_type != 'U:U':
                continue

            candidate_name = '_'.join([candidate_info.comp_name,
                                       candidate_info.candidate_id])

            if candidate_name in processed_candidates:
                continue

            try:
                process_candidate(candidate_info, options.db,
                                  options.prev_epoch, gap_stats_csv_b,
                                  processed_comps,
                                  gap_stats_csv_u)
            except Exception as e:
                print('Couldn\'t process candidate:', candidate_name,
                      'reason:', e, flush=True)
