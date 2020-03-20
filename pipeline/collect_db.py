import json
import math
import operator

import requests

from Bio.PDB import PDBList, PDBParser, PDBIO, Selection, Polypeptide, \
    PPBuilder
from pandas import read_csv
import os
import shutil
from xml.etree import ElementTree
from functools import reduce
from Bio import pairwise2
from collections import defaultdict
import sys
import signal
import functools

import alignments

complexes = []

PDB_ID = 'pdb'
H_CHAIN = 'Hchain'
L_CHAIN = 'Lchain'
ANTIGEN_CHAIN = 'antigen_chain'
ANTIGEN_TYPE = 'antigen_type'
ANTIGEN_HET_NAME = 'antigen_het_name'
STRUCTURE = 'structure'

NA = 'NA'

DB_PATH = 'data'
DOT_PDB = '.pdb'
DOT_FASTA = '.fasta'

MISMATCHED_LOG = 'mismatched.log'

AG = 'AG'
AB = 'AB'

CHAINS_SEPARATOR = '+'

peptides_builder = PPBuilder()


class HandlerError(RuntimeError):
    pass


def extract_seq(chain):
    seq = ''

    for x in peptides_builder.build_peptides(chain):
        seq += str(x.get_sequence())

    return seq


def pdb_seqs_for_structure(struct, chain_ids):
    pdb_seqs = []

    for chain_id in chain_ids:
        for model in struct:
            for chain in model:
                if chain.get_id() == chain_id:
                    pdb_seqs.append(extract_seq(chain))

    return pdb_seqs


def fetch_struct(pdb_id):
    curl = 'https://files.rcsb.org/download/{}.pdb'. \
        format(pdb_id)

    path_to_tmp = os.path.join(DB_PATH, pdb_id + DOT_PDB)

    if os.path.exists(path_to_tmp):
        return path_to_tmp

    r = get_while_true(curl)

    if r is None:
        return None

    with open(path_to_tmp, 'w') as f:
        f.write(r)

    return path_to_tmp


def with_timeout(timeout=None):
    def inner(f):
        def handler(*args):
            raise HandlerError()

        @functools.wraps(f)
        def inner_inner(*args, **kwargs):
            if timeout:
                signal.signal(signal.SIGALRM, handler)
                signal.alarm(timeout)

            try:
                res = f(*args, **kwargs)
            except Exception as e:
                return None

            if timeout:
                signal.alarm(0)

            return res

        return inner_inner

    return inner


memo_db = {}


def memoize(f):
    @functools.wraps(f)
    def inner(*args, **kwargs):
        new_args = frozenset(map(
            lambda x: x if not isinstance(x, dict) and not isinstance(x, list)
            else (json.dumps(x, sort_keys=True) if not isinstance(x, list)
                  else frozenset(x)), args))
        key = (f.__name__, new_args)

        if key in memo_db.keys():
            return memo_db[key]

        res = f(*args, **kwargs)
        memo_db[key] = res

        return res

    return inner


@with_timeout(timeout=100)
@memoize
def get_while_true(curl):
    not_finished = True

    content = None

    while not_finished:
        try:
            print('getting', curl, flush=True)
            res = requests.get(curl)
            content = res.content.decode('utf-8')

            if not content:
                continue

            if '404 Not Found' in content:
                return None

            if content.startswith('<!DOCTYPE'):
                continue

            not_finished = False
        except HandlerError:
            return None
        except Exception:
            pass

    return content


@with_timeout(timeout=100)
@memoize
def post_while_true(url, json):
    not_finished = True

    content = None

    while not_finished:
        try:
            print('posting', url, json, flush=True)
            res = requests.post(url, json)
            content = res.content.decode('utf-8')

            if not content:
                continue

            if '404 Not Found' in content:
                return None

            if content.startswith('<!DOCTYPE'):
                continue

            not_finished = False
        except HandlerError:
            return None
        except Exception:
            pass

    return content


def is_obsolete(pdb_id):
    curl = 'https://www.rcsb.org/pdb/rest/getEntityInfo?structureId={}' \
        .format(pdb_id)

    r = get_while_true(curl)
    xml = ElementTree.fromstring(r)

    for child in xml:
        if child.tag == 'obsolete':
            return True

    return False


def form_comp_name(pdb_id, ab_chains, ag_chains):
    ab_names = list(map(lambda x: x if x else '', ab_chains))
    comp_name = pdb_id + '_' + CHAINS_SEPARATOR.join(ab_names) + '|' + \
                CHAINS_SEPARATOR.join(ag_chains)
    return comp_name


def comp_name_to_pdb_and_chains(comp_name):
    [pdb_id, chains] = comp_name.split('_')
    ab_chains_s, ag_chains_s = chains.split('|')

    ab_chains = ab_chains_s.split(CHAINS_SEPARATOR)
    ag_chains = ag_chains_s.split(CHAINS_SEPARATOR)

    return pdb_id, ab_chains, ag_chains


def get_real_seqs(struct, chain_ids_to_seqs):
    chain_ids = []
    fasta_seqs = []

    for k, v in chain_ids_to_seqs:
        chain_ids.append(k)
        fasta_seqs.append(v)

    pdb_seqs = pdb_seqs_for_structure(struct, chain_ids)

    real_seqs = []

    for fasta_seq, pdb_seq in zip(fasta_seqs, pdb_seqs):
        alignment_n = alignments.align_possibly_gapped_sequence_on_its_complete_version(
            pdb_seq, fasta_seq)

        if not alignment_n:
            real_seqs.append(None)

        alignment = alignment_n[0]

        _, _, real_seq = cut_alignments(alignment[0],
                                        alignment[1])

        real_seqs.append(real_seq)

    return real_seqs


class Complex:
    pdb_parser = PDBParser()

    def __init__(self, pdb_id, h_chain, l_chain, antigen_chain,
                 antigen_het_name):
        self.pdb_id = pdb_id
        self.antibody_h_chain = h_chain
        self.antibody_l_chain = l_chain

        self.is_vhh = self.antibody_l_chain is None

        # if chain ids of antibody's chains are equal up to case,
        # it means that antibody has only one chain
        if self.antibody_h_chain and self.antibody_l_chain and \
                self.antibody_h_chain.upper() == self.antibody_l_chain.upper():
            self.antibody_h_chain = self.antibody_h_chain.upper()
            self.antibody_l_chain = None

        self.antigen_chains = antigen_chain
        self.antigen_het_name = antigen_het_name

        self.antibody_chains = [self.antibody_h_chain,
                                self.antibody_l_chain] if not self.is_vhh else [
            self.antibody_h_chain]

        self.struct = self.pdb_parser.get_structure(self.pdb_id,
                                                    fetch_struct(self.pdb_id))

        self.comp_name = form_comp_name(self.pdb_id, self.antibody_chains,
                                        self.antigen_chains)

        self.complex_dir_path = os.path.join(DB_PATH, self.pdb_id)

        self.antigen_seqs = get_real_seqs(self.struct,
                                          [(x, self._fetch_sequence(x)) for x
                                           in
                                           self.antigen_chains])

        self.antibody_h_seq = None

        if self.antibody_h_chain:
            self.antibody_h_seq = \
                get_real_seqs(self.struct, [
                    (self.antibody_h_chain,
                     self._fetch_sequence(self.antibody_h_chain))])[0]

        self.antibody_l_seq = None

        if self.antibody_l_chain:
            self.antibody_l_seq = \
                get_real_seqs(self.struct, [(self.antibody_l_chain,
                                             self._fetch_sequence(
                                                 self.antibody_l_chain))])[0]

        self.antibody_seqs = [self.antibody_h_seq,
                              self.antibody_l_seq] if not self.is_vhh else [
            self.antibody_h_seq]

    def has_unfetched_sequences(self):
        return list(filter(lambda x: x is None,
                           self.antigen_chains + self.antibody_chains))

    def _fetch_sequence(self, chain_id):
        fasta_path = os.path.join(self.complex_dir_path,
                                  self.pdb_id + '_' + chain_id + DOT_FASTA)

        print('fetching', fasta_path, flush=True)

        if os.path.exists(fasta_path):
            with open(fasta_path, 'r') as f:
                fasta = f.readlines()

            if len(fasta) < 2:
                os.remove(fasta_path)
                return self._fetch_sequence(chain_id)

            return fasta[1]

        fasta = ['> ' + self.pdb_id + ':' + chain_id,
                 fetch_sequence(self.pdb_id, chain_id)]

        if fasta[1] is None:
            return None

        if not os.path.exists(self.complex_dir_path):
            os.mkdir(self.complex_dir_path)

        with open(fasta_path, 'w') as f:
            f.write(fasta[0] + '\n' + fasta[1])

        return fasta[1]


def fetch_all_sequences(pdb_id):
    url = 'https://www.rcsb.org/pdb/download/downloadFastaFiles.do'
    r = post_while_true(url, {'structureIdList': pdb_id,
                              'compressionType': 'uncompressed'})

    seqs = []

    for line in r.split():
        if line.startswith('>'):
            seqs.append([line[6], ''])
        else:
            if not seqs:
                print('bad line:', line, 'in', r, flush=True)
                return fetch_all_sequences(pdb_id)

            seqs[-1][1] += line

    return {y[0]: y[1] for y in seqs}


def fetch_sequence(pdb_id, chain_id):
    seqs = fetch_all_sequences(pdb_id)

    if chain_id in seqs:
        return seqs[chain_id]

    return None


def sub_nan(val):
    if isinstance(val, float) and math.isnan(val):
        return None
    return val


def get_bound_complexes(sabdab_summary_df, to_accept=None, p=None):
    complexes = []

    obsolete = {}

    with open('obsolete.log', 'a+') as obsolete_log:
        obsolete_log.seek(0)

        for line in obsolete_log.readlines():
            key, value = line.strip().split(',')
            obsolete[key] = bool(int(value))

        counter = -1

        allowed_types_of_antigen = ['protein', 'protein | protein',
                                    'protein | protein | protein']

        for _, row in sabdab_summary_df.iterrows():
            counter += 1

            if p and not (p[0] <= counter < p[1]):
                continue

            if sub_nan(row[ANTIGEN_TYPE]) and row[ANTIGEN_TYPE] in \
                    allowed_types_of_antigen:
                if to_accept and row[PDB_ID].upper() not in to_accept:
                    continue

                is_vhh_l = sub_nan(row[H_CHAIN]) is None and sub_nan(
                    row[L_CHAIN]) is not None

                if row[PDB_ID] in obsolete.keys():
                    if obsolete[row[PDB_ID]]:
                        continue
                else:
                    is_obs = is_obsolete(row[PDB_ID])

                    obsolete_log.write(
                        '{},{}\n'.format(row[PDB_ID], int(is_obs)))
                    obsolete_log.flush()

                    if is_obs:
                        continue

                antigen_chains = row[ANTIGEN_CHAIN].split(' | ')

                if is_vhh_l:
                    new_complex = Complex(
                        row[PDB_ID], sub_nan(row[L_CHAIN]),
                        sub_nan(row[H_CHAIN]),
                        antigen_chains, sub_nan(row[ANTIGEN_HET_NAME]))
                else:
                    new_complex = Complex(
                        row[PDB_ID], sub_nan(row[H_CHAIN]),
                        sub_nan(row[L_CHAIN]),
                        antigen_chains, sub_nan(row[ANTIGEN_HET_NAME]))

                if new_complex.has_unfetched_sequences():
                    print('Has unfetched sequences:', row[PDB_ID])
                    continue

                complexes.append(new_complex)
            else:
                print('Not protein-protein complex:', row[PDB_ID])

    return complexes


class Candidate:
    def __init__(self, pdb_id, chain_ids):
        self.pdb_id = pdb_id
        self.chain_ids = chain_ids

    def __str__(self):
        return str((self.pdb_id, self.chain_ids))

    def __repr__(self):
        return str((self.pdb_id, self.chain_ids))


def cut_alignments(query_alignment, target_alignment):
    match_ids = []

    for i in range(len(query_alignment)):
        if query_alignment[i] != '-' and query_alignment[i] == \
                target_alignment[i]:
            match_ids.append(i)

    if not match_ids:
        return [], [], []

    first_match_id = match_ids[0]
    last_match_id = match_ids[-1]

    return match_ids, query_alignment[first_match_id: last_match_id + 1], \
           target_alignment[first_match_id: last_match_id + 1]


def calc_mismatches_stat(query_seq, target_seq):
    alignment_list = alignments.subsequence_without_gaps(query_seq, target_seq)

    if not alignment_list:
        return 0, len(query_seq), len(query_seq)

    alignment = alignment_list[0]

    mismatches_count = 0

    query_alignment = alignment[0]
    target_alignment = alignment[1]

    match_ids, query_alignment, target_alignment = cut_alignments(
        query_alignment, target_alignment)

    if not match_ids:
        return 0, len(query_seq), len(query_seq)

    cur_miss_len = 0
    max_miss_len = 0

    for i in range(len(query_alignment)):
        if query_alignment[i] != target_alignment[i]:
            cur_miss_len += 1
            mismatches_count += 1
        else:
            max_miss_len = max(max_miss_len, cur_miss_len)
            cur_miss_len = 0

    return len(match_ids), mismatches_count, max_miss_len


def is_subsequence_of(query_seq, target_seq, is_ab=True):
    min_intersection_len = int(0.9 * len(query_seq))
    max_miss_cutoff = max(10, int(0.03 * len(query_seq))) if is_ab else int(
        0.05 * len(query_seq))

    matches_count, mismatches_count, max_miss_len = calc_mismatches_stat(
        query_seq,
        target_seq)

    return matches_count + mismatches_count >= min_intersection_len and (
            not is_ab or max_miss_len < 3) and mismatches_count <= \
           max_miss_cutoff


# IDEA: maybe, make query seq the subseq of query seq that contains
#       all near-interface residues.
def is_match(query_seq, hit_alignment,
             is_ab=True):
    if query_seq == hit_alignment:
        return True

    hit_with_removed_gaps = hit_alignment.replace('-', '')

    if len(hit_with_removed_gaps) < int(0.9 * len(query_seq)):
        return False

    return is_subsequence_of(query_seq, hit_with_removed_gaps, is_ab=is_ab)


def get_blast_data(pdb_id, chain_id, seq, is_ab):
    curl = 'https://www.rcsb.org/pdb/rest/getBlastPDB2?structureId' \
           '={}&chainId={}&eCutOff=10.0&matrix=BLOSUM62&outputFormat=XML'. \
        format(pdb_id, chain_id)

    r = get_while_true(curl)
    xml = ElementTree.fromstring(r)

    res = []

    for child in xml:
        for iteration in child:
            for iteration_data in iteration:
                for hit in iteration_data:
                    if hit.tag != 'Hit':
                        continue

                    hit_def = hit.find('Hit_def')
                    hit_def_parts = hit_def.text.split('|')[0].split(':')

                    hit_pdb_id = hit_def_parts[0]

                    hit_chain_ids = [x for x in hit_def_parts[2].split(',')]

                    for hsp in hit.find('Hit_hsps'):
                        hsp_hseq = hsp.find('Hsp_hseq').text

                        good_chain_ids = []

                        for hit_chain_id in hit_chain_ids:
                            if is_match(seq, hsp_hseq,
                                        is_ab=is_ab):
                                good_chain_ids.append(hit_chain_id)

                        if good_chain_ids:
                            res.append(Candidate(hit_pdb_id, good_chain_ids))

    return res


def retrieve_uniprot_ids(pdb_id):
    url = 'https://www.uniprot.org/uploadlists/'
    r = post_while_true(url, {'from': 'PDB_ID',
                              'to': 'ACC',
                              'format': 'tab',
                              'query': pdb_id
                              })

    res = []

    for line in r.split('\n')[1:-1]:
        res.append(line.split('\t')[1])

    return res


def retrieve_names(pdb_id):
    curl = 'https://www.rcsb.org/pdb/rest/describeMol?structureId={}' \
        .format(pdb_id)

    r = get_while_true(curl)
    xml = ElementTree.fromstring(r)

    res = []

    for child in xml:
        for polymer in child:
            for attr in polymer:
                if attr.tag == 'polymerDescription':
                    res.append(attr.attrib['description'])

    return res


def retrieve_resolution(pdb_id):
    curl = 'https://www.rcsb.org/pdb/rest/getEntityInfo?structureId={}' \
        .format(pdb_id)

    r = get_while_true(curl)
    xml = ElementTree.fromstring(r)

    res = []

    for pdb in xml:
        try:
            res.append(pdb.attrib['resolution'])
        except Exception:
            # if there's no info about resolution,
            # then we consider it to be bad
            res.append(100)

    return float(res[0])


def check_names(names):
    if len(list(frozenset(names))) == 1:
        return True

    split_names = list(
        map(lambda x: list(map(lambda t: t.upper(), x.split())), names))

    if len(list(frozenset(map(lambda x: len(x), split_names)))) != 1:
        return False

    common_set = set(split_names[0])

    for x in split_names[1:]:
        common_set &= set(x)

    uncommon_set = set([])

    for x in split_names:
        for y in x:
            if y not in common_set:
                uncommon_set.add(y.upper())

    unknown_list = list(uncommon_set)

    for x in uncommon_set:
        if all(map(lambda t: t not in x, ['HEAVY', 'LIGHT'])):
            return False

    return len(unknown_list) == 2


def check_unbound(candidate_pdb_id, candidate_chain_ids, query_seqs, is_ab):
    candidate_seqs_dict = fetch_all_sequences(candidate_pdb_id)
    candidate_seqs = list(map(lambda x: candidate_seqs_dict[x],
                              candidate_chain_ids))

    for seq in candidate_seqs:
        if 'X' in seq:
            return None

    for i in range(len(query_seqs)):
        if not is_subsequence_of(query_seqs[i], candidate_seqs[i],
                                 is_ab=is_ab):
            return None

    c1 = check_names(retrieve_names(candidate_pdb_id))

    if not c1:
        return None

    return Candidate(candidate_pdb_id, candidate_chain_ids)


def sort_and_take_unbound(unbound_candidates):
    unbound_candidates.sort(key=lambda x: -int(x[0]))
    return unbound_candidates[:50]


def all_id_sets(ls, n):
    res = []

    for i in range(n):
        acc = []

        for x in ls:
            acc.append(x[i])

        res.append(acc)

    return res


def find_unbound_structure(pdb_id, chain_ids, seqs, is_ab):
    candidates = [get_blast_data(pdb_id, chain_id, seq, is_ab) for
                  chain_id, seq in
                  zip(chain_ids, seqs)]

    candidates_dicts = [{x.pdb_id: x.chain_ids for x in l} for l in candidates]

    pdb_ids_in_intersection_prep = reduce(operator.and_,
                                          [set([x.pdb_id for x in candidate])
                                           for
                                           candidate in candidates])

    unbound_candidates = \
        sort_and_take_unbound(list(pdb_ids_in_intersection_prep))

    res = []

    for candidate_id in unbound_candidates:
        if candidate_id.upper() == pdb_id.upper():
            continue

        candidate_chain_idss = []

        for x in candidates_dicts:
            candidate_chain_idss.append(x[candidate_id])

        # this works on a hunch
        all_sets_of_candidate_ids = all_id_sets(candidate_chain_idss,
                                                len(candidate_chain_idss[0]))

        for set_of_chain_ids in all_sets_of_candidate_ids:
            res_for_candidate = check_unbound(candidate_id, set_of_chain_ids,
                                              seqs, is_ab)

            if res_for_candidate:
                res.append(res_for_candidate)

    return res


def sort_and_take_ress(unbound_ress):
    unbound_ress.sort(key=lambda x: retrieve_resolution(x.pdb_id))

    taken_ids = set()
    res = []

    for candidate in unbound_ress:
        if candidate.pdb_id not in taken_ids:
            res.append(candidate)
            taken_ids.add(candidate.pdb_id)

    return res[:5]


def find_unbound_conformations(complex):
    unbound_antigen_valid_candidates = \
        find_unbound_structure(complex.pdb_id, complex.antigen_chains,
                               complex.antigen_seqs, False)

    print('unbound antigen:', unbound_antigen_valid_candidates, flush=True)

    unbound_antibody_valid_candidates = \
        find_unbound_structure(complex.pdb_id,
                               complex.antibody_chains, complex.antibody_seqs,
                               True)

    print('unbound antibody:', unbound_antibody_valid_candidates, flush=True)

    return sort_and_take_ress(unbound_antigen_valid_candidates), \
           sort_and_take_ress(unbound_antibody_valid_candidates)


structures_summary = read_csv('data/sabdab_summary_all.tsv',
                              sep='\t')

test_structures = [('1AHW', '1FGN', '1TFH'),
                   ('1BVK', '1BVL', '3LZT'),
                   ('1DQJ', '1DQQ', '3LZT'),
                   ('1E6J', '1E6O', '1A43'),
                   ('1JPS', '1JPT', '1TFH'),
                   ('1MLC', '1MLB', '3LZT'),
                   ('1VFB', '1VFA', '8LYZ'),
                   ('1WEJ', '1QBL', '1HRC'),
                   ('2FD6', '2FAT', '1YWH'),
                   ('2VIS', '1GIG', '2VIU'),
                   ('2VXT', '2VXU', '1J0S'),
                   ('2W9E', '2W9D', '1QM1'),
                   ('3EOA', '3EO9', '3F74'),
                   ('3HMX', '3HMW', '1F45'),
                   ('3MXW', '3MXV', '3M1N'),
                   ('3RVW', '3RVT', '3F5V'),
                   ('4DN4', '4DN3', '1DOL'),
                   ('4FQI', '4FQH', '2FK0'),
                   ('4G6J', '4G5Z', 'H5N1'),
                   ('4G6M', '4G6K', '4I1B'),
                   ('4GXU', '4GXV', '4I1B')]


def run_zlab_test():
    with open(MISMATCHED_LOG, 'w') as f:
        f.write(
            'bound_id,unbound_id,bound_chain,unbound_chain,mismatches_count,' +
            'len_diff\n')

    comps = get_bound_complexes(structures_summary,
                                list(map(lambda x: x[0], test_structures)))

    for pdb_id, unbound_antibody_id, unbound_antigen_id in test_structures:
        print('processing', pdb_id, flush=True)

        comps_found = list(filter(lambda x: x.pdb_id.upper() == pdb_id, comps))

        for comp in comps_found:
            unbound_antigen_candidates, unbound_antibody_candidates = \
                find_unbound_conformations(comp)

            print(comp.comp_name)

            print('antigen', 'expected:', unbound_antigen_id, 'got:',
                  unbound_antigen_candidates, flush=True)
            print('antibody', 'expected:', unbound_antibody_id, 'got:',
                  unbound_antibody_candidates)

            if unbound_antigen_id not in list(
                    map(lambda x: x.pdb_id, unbound_antigen_candidates)):
                print('MISMATCH! in antigen', flush=True)

            if unbound_antibody_id not in list(
                    map(lambda x: x.pdb_id, unbound_antibody_candidates)):
                print('MISMATCH! in antibody', flush=True)


@with_timeout(timeout=10)
def retrieve_pdb(pdb_id):
    url = 'https://files.rcsb.org/download/{}.pdb'.format(pdb_id)
    path_to_tmp = os.path.join(DB_PATH, pdb_id + DOT_PDB)

    if os.path.exists(path_to_tmp):
        return path_to_tmp

    res = post_while_true(url, {})

    with open(path_to_tmp, 'w') as f:
        f.write(res)

    return path_to_tmp


def remove_if_contains(path, s):
    for file in os.listdir(path):
        if s in file:
            os.remove(os.path.join(path, file))


def collect_unbound_structures(overwrite=True, p=None, to_accept=None):
    comps = get_bound_complexes(structures_summary, p=p,
                                to_accept=to_accept)

    processed = set()

    w_or_a = 'w' if overwrite else 'a'
    processed_open_mode = 'w' if overwrite else 'a+'

    with open('not_processed.log', w_or_a) as not_processed, open(
            'processed.log', processed_open_mode) as processed_log, \
            open('unbound_data.csv',
                 w_or_a) as unbound_data_csv:

        if overwrite:
            unbound_data_csv.write('pdb_id,comp_name,candidate_type,' +
                                   'candidate_pdb_id,candidate_chain_ids\n')
            unbound_data_csv.flush()

        if not overwrite:
            processed_log.seek(0)

            for processed_complex in processed_log.readlines():
                processed.add(processed_complex.strip())

        print('Complexes to process:', len(comps))

        for comp in comps:
            if comp.comp_name in processed:
                continue

            try:
                print('processing:', comp.comp_name, flush=True)

                unbound_antigen_candidates, unbound_antibody_candidates = \
                    find_unbound_conformations(comp)

                remove_if_contains(comp.complex_dir_path, AG)
                remove_if_contains(comp.complex_dir_path, AB)

                def helper_writer(candidates, suf):
                    counter = 0
                    for candidate in candidates:
                        print('candidate:', candidate, flush=True)

                        unbound_data_csv.write(
                            '{},{},{},{},{}\n'.format(comp.pdb_id,
                                                      comp.comp_name, suf,
                                                      candidate.pdb_id,
                                                      ':'.join(
                                                          candidate.
                                                              chain_ids)))
                        unbound_data_csv.flush()

                        counter += 1

                helper_writer(unbound_antigen_candidates, AG)
                helper_writer(unbound_antibody_candidates, AB)

                processed.add(comp.comp_name)
                processed_log.write(comp.comp_name + '\n')
                processed_log.flush()
            except Exception as e:
                not_processed.write('{}: {}\n'.format(comp.comp_name, e))
                not_processed.flush()


if __name__ == '__main__':
    if sys.argv and sys.argv[0] == 'test':
        run_zlab_test()
    else:
        p = list(filter(lambda x: x.startswith('--range='), sys.argv))

        if p:
            rest = p[0][8:].strip('(').strip(')').split(',')
            p = (int(rest[0]), int(rest[1]))
        else:
            p = None

        to_accept = None

        if os.path.exists('uus.txt'):
            with open('uus.txt', 'r') as f:
                to_accept = list(map(lambda x: x[:4].upper(), f.readlines()))

        collect_unbound_structures(overwrite=len(
            list(filter(lambda x: x == 'continue', sys.argv))) == 0, p=p,
                                   to_accept=to_accept)

# 1jps,1jps_H+L|T,AB,1JPT,H:L
# BAD: '4I2X', '5NH3', '1S78', '4XT1', '6OGE' (clone), '5WT9'
# GOOD: '6RCV', '6RCU',
