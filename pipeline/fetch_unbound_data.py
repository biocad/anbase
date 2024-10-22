import json
import math
import operator
import random
import time
import traceback
from multiprocessing import Process, Queue
from multiprocessing.pool import Pool
import re

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
from urllib.parse import quote_plus
import json

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

            if '429 Too Many Requests' in content:
                roll = random.random()
                time.sleep(roll)

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

            if '429 Too Many Requests' in content:
                roll = random.random()
                time.sleep(roll)

            if content.startswith('<!DOCTYPE'):
                continue

            not_finished = False
        except HandlerError:
            return None
        except Exception:
            pass

    return content


def is_obsolete(pdb_id):
    # curl = 'https://www.rcsb.org/pdb/rest/getEntityInfo?structureId={}' \
    #     .format(pdb_id)
    # 
    # r = get_while_true(curl)
    # xml = ElementTree.fromstring(r)
    # 
    # for child in xml:
    #     if child.tag == 'obsolete':
    #         return True

    return False


def form_comp_name(pdb_id, ab_chains, ag_chains):
    ab_names = list(map(lambda x: x if x else '', ab_chains))
    comp_name = pdb_id + '_' + CHAINS_SEPARATOR.join(ab_names) + '-' + \
                CHAINS_SEPARATOR.join(ag_chains)
    return comp_name


def comp_name_to_pdb_and_chains(comp_name):
    [pdb_id, chains] = comp_name.split('_')
    ab_chains_s, ag_chains_s = chains.split('-')

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
            continue

        alignment = alignment_n[0]

        _, _, real_seq = cut_alignments(alignment[0], alignment[1])

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
            chain_ids_to_seqs = [(self.antibody_l_chain, self._fetch_sequence(self.antibody_l_chain))]
            real_seqs = get_real_seqs(self.struct, chain_ids_to_seqs)
            self.antibody_l_seq = real_seqs[0]

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

def parse_freaking_chain_names(chains_line):
    chains_line_tmp = chains_line

    chain_names_auth = [m.replace('auth ', '') for m in re.findall('auth .', chains_line_tmp)]

    chains_line_tmp = re.sub('\[auth..\]', '', chains_line_tmp) # remove all strings like '[auth A]'
    chains_line_tmp = chains_line_tmp.replace('Chains','').replace('Chain','') # remove all 'Chains' and 'Chain'
    chains_line_tmp = chains_line_tmp.replace(' ', '') # remove spaces and 
    chains_names = chains_line_tmp.split(',') # just split by ','

    all_chain_names = chain_names_auth + chains_names

    return all_chain_names

def fetch_all_sequences(pdb_id, mol_names_res=None):
    url = f'https://www.rcsb.org/fasta/entry/{pdb_id}'
    r = get_while_true(url)

    seqs = {}
    chain_names = []

    for line in r.split('\n'):
        if line.startswith('>'):
            [_, chains, mol_name] = line.strip().split('|')[:3]

            if mol_names_res is not None:
                mol_names_res.append(mol_name)

            chain_names = parse_freaking_chain_names(chains)

            for chain_name in chain_names:
                seqs[chain_name] = ''
        else:
            if not chain_names:
                print('bad line:', line, 'in', r, flush=True)
                return fetch_all_sequences(pdb_id)

            for chain_name in chain_names:
                seqs[chain_name] += line

    return seqs

def fetch_all_sequences_for_entity(pdb_id, entity_id):
    url = f'https://www.rcsb.org/fasta/entry/{pdb_id}'
    r = get_while_true(url)

    seqs = {}
    chain_names = []

    filling = False

    for line in r.split('\n'):
        if line.startswith('>'):
            [pdb_entity, chains] = line.strip().split('|')[:2]

            this_entity_id = pdb_entity.split('_')[1]

            if this_entity_id != entity_id:
                filling = False
                continue

            chain_names = parse_freaking_chain_names(chains)

            for chain_name in chain_names:
                seqs[chain_name] = ''

            filling = True
        else:
            if not filling:
                continue

            if not chain_names:
                print('bad line:', line, 'in', r, flush=True)
                return fetch_all_sequences(pdb_id)

            for chain_name in chain_names:
                seqs[chain_name] += line

    return seqs


def fetch_sequence(pdb_id, chain_id):
    seqs = fetch_all_sequences(pdb_id)

    if chain_id in seqs:
        return seqs[chain_id]

    return None


def sub_nan(val):
    if isinstance(val, float) and math.isnan(val):
        return None
    return val


def get_bound_complexes(run_id, sabdab_summary_df, to_accept=None, p=None):
    complexes = []

    obsolete = {}

    with open('obsolete_{}.log'.format(run_id), 'a+') as obsolete_log:
        obsolete_log.seek(0)

        for line in obsolete_log.readlines():
            key, value = line.strip().split(',')
            obsolete[key] = bool(int(value))

        counter = -1

        allowed_types_of_antigen = ['protein', 'protein | protein',
                                    'protein | protein | protein']

        for _, row in sabdab_summary_df.iterrows():
            try:
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
            except Exception as e:
                print(f'WARNING: Complex {row[PDB_ID]} is not read:', e)
                traceback.print_tb(e.__traceback__)
                # raise e

    return complexes


class Candidate:
    def __init__(self, pdb_id, chain_ids):
        self.pdb_id = pdb_id
        self.chain_ids = chain_ids

    def __str__(self):
        return str((self.pdb_id, self.chain_ids))

    def __repr__(self):
        return str((self.pdb_id, self.chain_ids))

def is_subsequence_of(query_seq, target_seq):
  identity = alignments.calc_identity(query_seq, target_seq)

  return identity == 1.0

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


def calc_mismatches(query_seq, target_seq):
    alignment_list = alignments.subsequence_without_gaps(query_seq, target_seq)

    if not alignment_list:
        return len(query_seq)

    alignment = alignment_list[0]

    mismatches_count = 0

    query_alignment = alignment[0]
    target_alignment = alignment[1]

    match_ids, query_alignment, target_alignment = cut_alignments(
        query_alignment, target_alignment)

    if not match_ids:
        return len(query_seq)

    for i in range(len(query_alignment)):
      if query_alignment[i] != target_alignment[i]:
        mismatches_count += 1

    return mismatches_count

class Hit:
  def __init__(self, hit):
    [self.pdb_id, entity_id] = hit['identifier'].split('_')
    self.chain_ids_to_seqs = fetch_all_sequences_for_entity(self.pdb_id,
                                                                entity_id)

def process_hit(seq, hit):
    res = []
    good_chain_ids = []

    for chain_id, hit_seq in hit.chain_ids_to_seqs.items():

        if is_subsequence_of(hit_seq, seq):
            good_chain_ids.append(chain_id)

    if good_chain_ids:
        res.append(Candidate(hit.pdb_id, good_chain_ids))

    return res


def get_blast_data(pdb_id, chain_id, seq, is_ab):
    params = \
        {
            'query': {
                'type': 'terminal',
                'service': 'sequence',
                'parameters': {
                    'evalue_cutoff': 10,
                    'identity_cutoff': 0.9,
                    'target': 'pdb_protein_sequence',
                    'value': seq
                }
            },
            'request_options': {
                'scoring_strategy': 'sequence',
                "pager": {
                    "start": 0,
                    "rows": 100
                }

            },
            'return_type': 'polymer_entity'
        }

    curl = f'https://search.rcsb.org/rcsbsearch/v1/query?' \
           f'json={quote_plus(json.dumps(params, separators=(",", ":")))}'

    r = get_while_true(curl)

    res = []

    hits = []

    for hit in json.loads(r)['result_set']:
        hits.append((seq, Hit(hit)))

    # if NUMBER_OF_PROCESSES > 1:
    with Pool(NUMBER_OF_PROCESSES) as pool:
      for x in pool.starmap(process_hit, hits):
        res += x
    # else:
    #   for hit in hits:
    #     res += process_hit(hit)

    print(res)

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
    res = []

    fetch_all_sequences(pdb_id, mol_names_res=res)

    return res


def retrieve_resolution(pdb_id):
    curl = f'https://data.rcsb.org/rest/v1/core/entry/{pdb_id}'

    r = get_while_true(curl)

    unknown_method = 'UNKNOWN'

    info = json.loads(r)

    try:
        info = json.loads(r)

        resolution = int(info['pdbx_vrpt_summary']['pdbresolution'])
        method = info['rcsb_entry_info']['experimental_method']

        return resolution, method
    except Exception:
        # if there's no info about resolution,
        # then we consider it to be bad
        return 100, unknown_method


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
    unbound_ress.sort(key=lambda x: retrieve_resolution(x.pdb_id)[0])

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


def run_zlab_test(structures_summary):
    comps = get_bound_complexes('-1', structures_summary,
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


def collect_unbound_structures(run_id, structures_summary, overwrite=True, p=None, to_accept=None):
    comps = get_bound_complexes(run_id, structures_summary, p=p,
                                to_accept=to_accept)

    processed = set()

    w_or_a = 'w' if overwrite else 'a'
    processed_open_mode = 'w' if overwrite else 'a+'

    with open('not_processed_{}.log'.format(run_id),
              processed_open_mode) as not_processed, open(
        'processed_{}.log'.format(run_id),
        processed_open_mode) as processed_log, \
            open('unbound_data_{}.csv'.format(run_id),
                 w_or_a) as unbound_data_csv:

        if overwrite:
            unbound_data_csv.write('pdb_id,comp_name,candidate_type,' +
                                   'candidate_pdb_id,candidate_chain_ids\n')
            unbound_data_csv.flush()

        if not overwrite:
            processed_log.seek(0)

            for processed_complex in processed_log.readlines():
                processed.add(processed_complex.strip())

            not_processed.seek(0)

            for not_processed_complex in not_processed.readlines():
                processed.add(not_processed_complex.strip())

        print('Complexes to process:', len(comps))

        for comp in comps:
            if comp.comp_name in processed:
                continue

            try:
                print('processing:', comp.comp_name, flush=True)

                unbound_antigen_candidates, unbound_antibody_candidates = \
                    find_unbound_conformations(comp)

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
                print(e)
                traceback.print_tb(e.__traceback__)
                not_processed.write('{}: {}\n'.format(comp.comp_name, e))
                not_processed.flush()

if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('--sabdab-summary', default='sabdab_summary_all.tsv',
                      dest='sabdab_summary_file_path', metavar='SABDAB_SUMMARY',
                      help='Path to sabdab summary file')
    parser.add_option('--run-id', default='0',
                      dest='run_id',
                      metavar='RUN_ID',
                      help='ID of the current run [default: {}]'.
                      format('0'))
    parser.add_option('--is-test', default=False, dest='is_test',
                      metavar='IS_TEST',
                      help='Run in test mode [default: {}]'.
                      format('False'))
    parser.add_option('--range', default=None,
                      dest='range', metavar='RANGE',
                      help='Range of complexes to process from sabdab_'
                           'summary_all.tsv. [default: {}]'.format('None'))
    parser.add_option('--continue', default=False,
                      dest='cont', metavar='CONTINUE',
                      help='Whether to continue execution of a script from '
                           'the cached place [default: {}]'.format('False'))
    parser.add_option('--number-of-processes', default=3,
                      dest='n', metavar='N',
                      help='Number of paralles prcoesses')
    options, _ = parser.parse_args()

    structures_summary = read_csv(options.sabdab_summary_file_path,
                              sep='\t')

    NUMBER_OF_PROCESSES = int(options.n)

    if options.cont == 'True':
        cont = True
    else:
        cont = False

    if options.is_test:
        run_zlab_test(structures_summary)

    if options.range:
        [l, r] = options.range.replace(' ', '').strip('(').strip(')').split(
            ',')
        p = (int(l), int(r))
    else:
        p = None

    collect_unbound_structures(options.run_id, structures_summary, overwrite=not cont, p=p)
