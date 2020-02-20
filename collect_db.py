import json
import math
import operator

import requests

from Bio.PDB import PDBList, PDBParser, PDBIO, Selection, Polypeptide
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


class HandlerError(RuntimeError):
    pass


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


class Complex:
    pdb_parser = PDBParser()

    def __init__(self, pdb_id, h_chain, l_chain, antigen_chain,
                 antigen_het_name):
        self.pdb_id = pdb_id
        self.antibody_h_chain = h_chain
        self.antibody_l_chain = l_chain

        # if chain ids of antibody's chains are equal up to case,
        # it means that antibody has only one chain
        if self.antibody_h_chain and self.antibody_l_chain and \
                self.antibody_h_chain.upper() == self.antibody_l_chain.upper():
            self.antibody_h_chain = self.antibody_h_chain.upper()
            self.antibody_l_chain = None

        self.antigen_chains = antigen_chain
        self.antigen_het_name = antigen_het_name
        self.structure = None

        h_name = self.antibody_h_chain if self.antibody_h_chain else ''
        l_name = self.antibody_l_chain if self.antibody_l_chain else ''
        self.db_name = self.pdb_id + '_' + ''.join(
            [h_name, l_name] + self.antigen_chains)

        self.complex_dir_path = os.path.join(DB_PATH, self.pdb_id)

        self.antigen_seqs = [self._fetch_sequence(x) for x in
                             self.antigen_chains]

        self.antibody_h_seq = None

        if self.antibody_h_chain:
            self.antibody_h_seq = self._fetch_sequence(self.antibody_h_chain)

        self.antibody_l_seq = None

        if self.antibody_l_chain:
            self.antibody_l_seq = self._fetch_sequence(self.antibody_l_chain)

    def has_unfetched_sequences(self):
        return list(filter(lambda x: x is None, self.antigen_chains +
                           [self.antibody_h_chain, self.antibody_l_chain]))

    def load_structure(self):
        self.load_structure_from(os.path.join(self.complex_dir_path,
                                              self.db_name + DOT_PDB))

    def load_structure_from(self, path):
        try:
            self.structure = self.pdb_parser.get_structure(self.pdb_id, path)
        except Exception as e:
            print('Not loaded:', self.db_name, e)
            self.structure = None

    def is_bad_structure(self):
        return self.structure is None

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

    return list(map(lambda y: (y[0], y[1]), seqs))


def fetch_sequence(pdb_id, chain_id):
    seqs = fetch_all_sequences(pdb_id)

    potential_res = list(filter(lambda x: x[0] == chain_id, seqs))

    if potential_res:
        return potential_res[0][1]

    return None


def get_bound_complexes(sabdab_summary_df, to_accept=None, p=None):
    def sub_nan(val):
        if isinstance(val, float) and math.isnan(val):
            return None
        return val

    complexes = []

    obsolete = {}

    with open('obsolete.log', 'a+') as obsolete_log:
        obsolete_log.seek(0)

        for line in obsolete_log.readlines():
            key, value = line.strip().split(',')
            obsolete[key] = bool(int(value))

        counter = -1

        allowed_types_of_antigen = ['protein', 'peptide', 'protein | protein',
                                    'protein | protein | protein',
                                    'protein | peptide',
                                    'peptide | protein',
                                    'peptide | peptide | peptide',
                                    'protein | protein | peptide']

        for _, row in sabdab_summary_df.iterrows():
            counter += 1

            if p and not (p[0] <= counter < p[1]):
                continue

            if sub_nan(row[ANTIGEN_TYPE]) and row[ANTIGEN_TYPE] in \
                    allowed_types_of_antigen:
                if to_accept and row[PDB_ID].upper() not in to_accept:
                    continue

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
                new_complex = Complex(
                    row[PDB_ID], sub_nan(row[H_CHAIN]), sub_nan(row[L_CHAIN]),
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


def load_bound_complexes(complexes, load_structures=False):
    tmp_paths = set([])

    io = PDBIO()

    for comp in complexes:
        pdb_path = os.path.join(comp.complex_dir_path,
                                comp.db_name + DOT_PDB)

        if os.path.exists(pdb_path):
            if load_structures:
                comp.load_structure_from(pdb_path)
            print(comp.pdb_id, 'loaded', flush=True)
            continue

        # if os.path.exists(comp.complex_dir_path):
        #     shutil.rmtree(comp.complex_dir_path)

        if not os.path.exists(comp.complex_dir_path):
            os.mkdir(comp.complex_dir_path)

        tmp_path = retrieve_pdb(comp.pdb_id)

        if not tmp_path:
            continue

        comp.load_structure_from(tmp_path)

        if comp.is_bad_structure():
            continue

        needed_chain_ids = [x for x in [comp.antibody_h_chain,
                                        comp.antibody_l_chain] +
                            comp.antigen_chains if x]

        for model in comp.structure:
            for chain in model:
                if chain.get_id() not in needed_chain_ids:
                    model.detach_child(chain.get_id())

        io.set_structure(comp.structure)
        io.save(pdb_path)

        tmp_paths.add(tmp_path)

        print(comp.pdb_id, 'loaded', flush=True)

    for tmp_path in tmp_paths:
        os.remove(tmp_path)


def align_and_check(query_seq, target_seq, pdb_ids, chain_ids, write_log,
                    len_diff, is_ab=True):
    cut_off = int(0.03 * len(target_seq))

    alignment_list = pairwise2.align.localxs(query_seq, target_seq, -10, -10,
                                             penalize_end_gaps=False,
                                             one_alignment_only=True)

    if not alignment_list:
        return False

    alignment = alignment_list[0]

    mismatches_count = 0

    query_alignment = alignment[0]
    target_alignment = alignment[1]

    cur_miss_len = 0
    max_miss_len = 0

    for i in range(len(query_alignment)):
        if query_alignment[i] == '-' or target_alignment[i] == '-':
            max_miss_len = max(max_miss_len, cur_miss_len)
            cur_miss_len = 0
            continue

        if query_alignment[i] != target_alignment[i]:
            cur_miss_len += 1
            mismatches_count += 1
        else:
            max_miss_len = max(max_miss_len, cur_miss_len)
            cur_miss_len = 0

    if write_log and 0 < mismatches_count <= cut_off:
        with open(MISMATCHED_LOG, 'a') as f:
            f.write(','.join(
                [pdb_ids[0], pdb_ids[1], chain_ids[0].upper(),
                 chain_ids[1].upper(),
                 str(mismatches_count), str(len_diff)]) + '\n')

    return (not is_ab or max_miss_len < 3) and mismatches_count <= cut_off


def compare_query_and_hit_seqs(query_seq, hit_seq, pdb_ids, chain_ids,
                               write_log=False, is_ab=True):
    cut_off_half = int(0.1 * len(query_seq) / 2)
    len_diff = abs(len(query_seq) - len(hit_seq))

    c1 = align_and_check(query_seq, hit_seq, pdb_ids, chain_ids, write_log,
                         len_diff, is_ab=is_ab)

    if c1:
        return True

    c2 = align_and_check(query_seq[cut_off_half:-cut_off_half], hit_seq,
                         pdb_ids, chain_ids, write_log, len_diff, is_ab=is_ab)

    if c2:
        return True

    c3 = align_and_check(hit_seq[cut_off_half:-cut_off_half], query_seq,
                         pdb_ids, chain_ids, write_log, len_diff, is_ab=is_ab)

    if c3:
        return True

    c4 = align_and_check(query_seq[:-2 * cut_off_half], hit_seq, pdb_ids,
                         chain_ids, write_log, len_diff, is_ab=is_ab)

    if c4:
        return True

    c5 = align_and_check(hit_seq[:-2 * cut_off_half], query_seq, pdb_ids,
                         chain_ids, write_log, len_diff, is_ab=is_ab)

    if c5:
        return True

    c6 = align_and_check(query_seq[2 * cut_off_half:], hit_seq, pdb_ids,
                         chain_ids, write_log, len_diff, is_ab=is_ab)

    if c6:
        return True

    c7 = align_and_check(hit_seq[2 * cut_off_half:], query_seq, pdb_ids,
                         chain_ids, write_log, len_diff, is_ab=is_ab)

    return c7


def is_match(query_seq, query_alignment, hit_alignment, pdb_ids, chain_ids,
             is_ab=True):
    if query_seq == hit_alignment:
        return True

    query_with_stripped_gaps = query_alignment.strip('-')

    if '-' in query_with_stripped_gaps:
        return False

    hit_with_stripped_gaps = hit_alignment.strip('-')

    if '-' in hit_with_stripped_gaps:
        return False

    return compare_query_and_hit_seqs(query_seq, hit_with_stripped_gaps,
                                      pdb_ids, chain_ids, is_ab=is_ab)


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
                        hsp_qseq = hsp.find('Hsp_qseq').text
                        hsp_hseq = hsp.find('Hsp_hseq').text

                        if not is_match(seq, hsp_qseq, hsp_hseq,
                                        (pdb_id, hit_pdb_id),
                                        (chain_id, hit_chain_ids[0]),
                                        is_ab=is_ab):
                            continue

                        res.append(Candidate(hit_pdb_id, hit_chain_ids))

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


def check_unbound(pdb_id, chain_ids_and_seqs, query_pdb_id, is_ab):
    all_seqs_in_pdb = fetch_all_sequences(pdb_id)

    for _, seq in all_seqs_in_pdb:
        if 'X' in seq:
            return []

    chain_matches = defaultdict(list)

    for chain_id, chain_seq in chain_ids_and_seqs:
        for target_chain_id, seq in all_seqs_in_pdb:
            if compare_query_and_hit_seqs(chain_seq, seq,
                                          (query_pdb_id, pdb_id),
                                          (chain_id, target_chain_id),
                                          write_log=True, is_ab=is_ab):
                chain_matches[chain_id].append(target_chain_id)

    # we check that for every queried chain there is a matching chain in the
    # given pdb
    # also if names of all structures in pdb are different in no more than
    # one word (for example, 'my ab heavy chain' and 'my ab light chain)
    # it usually means that structures form one macromolecule,
    # hence their complex is unbound
    c1 = all(map(lambda x: len(chain_matches[x]) > 0, chain_matches.keys()))
    c2 = check_names(retrieve_names(pdb_id))

    if c1 and c2:
        res = []
        for i in range(len(chain_matches[chain_ids_and_seqs[0][0]])):
            res.append(Candidate(pdb_id, list(
                map(lambda x: chain_matches[x][i], chain_matches.keys()))))
        return res

    return []


def sort_and_take_unbound(unbound_candidates):
    unbound_candidates.sort(key=lambda x: -int(x[0]))
    return unbound_candidates[:50]


def find_unbound_structure(pdb_id, chain_ids, seqs, is_ab):
    candidates = [get_blast_data(pdb_id, chain_id, seq, is_ab) for
                  chain_id, seq in
                  zip(chain_ids, seqs)]

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

        res += check_unbound(candidate_id, list(zip(chain_ids, seqs)), pdb_id,
                             is_ab)

    return res


def sort_and_take_ress(unbound_ress):
    unbound_ress.sort(key=lambda x: retrieve_resolution(x.pdb_id))
    return unbound_ress[:5]


def find_unbound_conformations(complex):
    unbound_antigen_valid_candidates = \
        find_unbound_structure(complex.pdb_id, complex.antigen_chains,
                               complex.antigen_seqs, False)

    print('unbound antigen:', unbound_antigen_valid_candidates, flush=True)

    unbound_antibody_valid_candidates = \
        find_unbound_structure(complex.pdb_id,
                               [complex.antibody_h_chain,
                                complex.antibody_l_chain],
                               [complex.antibody_h_seq,
                                complex.antibody_l_seq], True)

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
    load_bound_complexes(comps)

    for pdb_id, unbound_antibody_id, unbound_antigen_id in test_structures:
        print('processing', pdb_id, flush=True)

        comps_found = list(filter(lambda x: x.pdb_id.upper() == pdb_id, comps))

        for comp in comps_found:
            comp.load_structure()

            unbound_antigen_candidates, unbound_antibody_candidates = \
                find_unbound_conformations(comp)

            print(comp.db_name)

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


def collect_unbound_structures(overwrite=True, p=None):
    comps = get_bound_complexes(structures_summary, p=p)
    load_bound_complexes(comps)

    tmp_paths = set()

    pdb_parser = PDBParser()
    io = PDBIO()

    processed = set()

    w_or_a = 'w' if overwrite else 'a'
    processed_open_mode = 'w' if overwrite else 'a+'

    with open('not_processed.log', w_or_a) as not_processed, open(
            'processed.log', processed_open_mode) as processed_log, \
            open('unbound_data.csv',
                 w_or_a) as unbound_data_csv:

        if overwrite:
            unbound_data_csv.write('pdb_id,db_name,type,candidate,' +
                                   'candidate_pdb_id,candidate_chain_names\n')
            unbound_data_csv.flush()

        if not overwrite:
            processed_log.seek(0)

            for processed_complex in processed_log.readlines():
                processed.add(processed_complex.strip())

        print('Complexes to process:', len(comps))

        for comp in comps:
            if comp.db_name in processed:
                continue

            try:
                comp.load_structure()

                if comp.is_bad_structure():
                    raise RuntimeError('Couldn\'t parse: ' + comp.db_name)

                print('processing:', comp.db_name, flush=True)

                unbound_antigen_candidates, unbound_antibody_candidates = \
                    find_unbound_conformations(comp)

                remove_if_contains(comp.complex_dir_path, AG)
                remove_if_contains(comp.complex_dir_path, AB)

                def helper_writer(candidates, suf):
                    counter = 0
                    for candidate in candidates:
                        print('candidate:', candidate, flush=True)

                        candidate_name = comp.pdb_id + '_' + \
                                         suf + '_u_' + str(counter)

                        path_to_candidate_pdb = os.path.join(
                            comp.complex_dir_path,
                            candidate_name + DOT_PDB)

                        tmp_path = retrieve_pdb(candidate.pdb_id)

                        if not tmp_path:
                            continue

                        structure = pdb_parser.get_structure(candidate.pdb_id,
                                                             tmp_path)

                        for model in comp.structure:
                            for chain in model:
                                if chain.get_id() not in candidate.chain_ids:
                                    model.detach_child(chain.get_id())

                        io.set_structure(structure)
                        io.save(path_to_candidate_pdb)

                        tmp_paths.add(tmp_path)

                        unbound_data_csv.write(
                            '{},{},{},{},{},{}\n'.format(comp.pdb_id,
                                                         comp.db_name, suf,
                                                         candidate_name,
                                                         candidate.pdb_id,
                                                         ':'.join(
                                                             candidate.
                                                                 chain_ids)))
                        unbound_data_csv.flush()

                        counter += 1

                helper_writer(unbound_antigen_candidates, AG)
                helper_writer(unbound_antibody_candidates, AB)

                processed.add(comp.db_name)
                processed_log.write(comp.db_name + '\n')
                processed_log.flush()
            except Exception as e:
                not_processed.write('{}: {}\n'.format(comp.db_name, e))
                not_processed.flush()

    for tmp_path in tmp_paths:
        os.remove(tmp_path)


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

        print(p)

        collect_unbound_structures(overwrite=len(
            list(filter(lambda x: x == 'continue', sys.argv))) == 0, p=p)
