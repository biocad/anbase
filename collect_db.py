import math
import operator

import requests

from Bio.PDB import PDBList, PDBParser, PDBIO, Selection, Polypeptide
from pandas import read_csv
import os
import shutil
from xml.etree import ElementTree
from functools import reduce

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


def get_while_true(curl):
    not_finished = True

    res = None

    while not_finished:
        try:
            res = requests.get(curl)
            not_finished = False
        except Exception:
            pass

    return res.content.decode('utf-8')


def post_while_true(url, json):
    not_finished = True

    res = None

    while not_finished:
        try:
            res = requests.post(url, json)
            not_finished &= False
        except Exception:
            pass

    return res.content.decode('utf-8')


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

        self.complex_dir_path = os.path.join(DB_PATH, self.pdb_id)

        self.antigen_seqs = [self._fetch_sequence(x) for x in
                             self.antigen_chains]

        self.antibody_h_seq = None

        if self.antibody_h_chain:
            self.antibody_h_seq = self._fetch_sequence(self.antibody_h_chain)

        self.antibody_l_seq = None

        if self.antibody_l_chain:
            self.antibody_l_seq = self._fetch_sequence(self.antibody_l_chain)

    def load_structure(self):
        self.load_structure_from(os.path.join(self.complex_dir_path,
                                              self.pdb_id + DOT_PDB))

    def load_structure_from(self, path):
        self.structure = self.pdb_parser.get_structure(self.pdb_id, path)

    def _fetch_sequence(self, chain_id):
        fasta_path = os.path.join(self.complex_dir_path,
                                  self.pdb_id + '_' + chain_id + DOT_FASTA)

        print('fetching', fasta_path)

        if os.path.exists(fasta_path):
            with open(fasta_path, 'r') as f:
                fasta = f.readlines()

            if len(fasta) < 2:
                os.remove(fasta_path)
                return self._fetch_sequence(chain_id)

            return fasta[1]

        fasta = ['> ' + self.pdb_id + ':' + chain_id,
                 fetch_sequence(self.pdb_id, chain_id)]

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
            seqs[-1][1] += line

    return list(map(lambda y: (y[0], y[1]), seqs))


def fetch_sequence(pdb_id, chain_id):
    seqs = fetch_all_sequences(pdb_id)

    return next(filter(lambda x: x[0] == chain_id, seqs))[1]


def get_bound_complexes(sabdab_summary_df, to_accept=None):
    def sub_nan(val):
        if isinstance(val, float) and math.isnan(val):
            return None
        return val

    complexes = []

    for _, row in sabdab_summary_df.iterrows():
        # if antigen's type is in lower case, it means that antigen is no good
        # for us, because it's a small molecule
        if sub_nan(row[ANTIGEN_TYPE]) and row[ANTIGEN_TYPE].islower():
            if to_accept and row[PDB_ID].upper() not in to_accept:
                continue

            antigen_chains = row[ANTIGEN_CHAIN].split(' | ')
            complexes.append(Complex(
                row[PDB_ID], sub_nan(row[H_CHAIN]), sub_nan(row[L_CHAIN]),
                antigen_chains,
                sub_nan(row[ANTIGEN_HET_NAME])))

    return complexes


class BLASTData:
    def __init__(self, pdb_id, chain_id):
        self.pdb_id = pdb_id
        self.chain_id = chain_id

    def __str__(self):
        return str((self.pdb_id, self.chain_id))


def load_bound_complexes(complexes, load_structures=False):
    with open('could_not_fetch.log', 'w') as could_not_fetch_log:
        pdb_list = PDBList()

        io = PDBIO()

        for comp in complexes:
            pdb_path = os.path.join(comp.complex_dir_path,
                                    comp.pdb_id + DOT_PDB)

            if os.path.exists(pdb_path):
                if load_structures:
                    comp.load_structure_from(pdb_path)
                print(comp.pdb_id, 'loaded')
                continue

            if os.path.exists(comp.complex_dir_path):
                shutil.rmtree(comp.complex_dir_path)

            os.mkdir(comp.complex_dir_path)

            ent_path = pdb_list.retrieve_pdb_file(comp.pdb_id,
                                                  file_format='pdb',
                                                  pdir=DB_PATH)

            if not os.path.exists(ent_path):
                print('Not written:', comp.pdb_id)
                print(comp.pdb_id, flush=True, file=could_not_fetch_log)
                continue

            comp.load_structure_from(ent_path)

            needed_chain_ids = [x for x in [comp.antibody_h_chain, comp.antibody_l_chain] +
                                comp.antigen_chains if x]

            for model in comp.structure:
                for chain in model:
                    if chain.get_id() not in needed_chain_ids:
                        model.detach_child(chain.get_id())

            io.set_structure(comp.structure)
            io.save(pdb_path)

            os.remove(ent_path)

            print(comp.pdb_id, 'loaded')


def compare_query_and_hit_seqs(query_seq, hit_seq):
    cut_off_half = int(0.05 * len(query_seq) / 2)

    c1 = query_seq in hit_seq
    c2 = hit_seq in query_seq and abs(len(hit_seq) - len(query_seq)) \
        <= cut_off_half * 2
    c3 = query_seq[cut_off_half:-cut_off_half] in hit_seq
    c4 = hit_seq[cut_off_half:-cut_off_half] in query_seq
    c5 = query_seq[:-2 * cut_off_half] in hit_seq
    c6 = hit_seq[:-2 * cut_off_half] in query_seq
    c7 = query_seq[2 * cut_off_half:] in hit_seq
    c8 = hit_seq[2 * cut_off_half:] in query_seq

    return c1 or c2 or c3 or c4 or c5 or c6 or c7 or c8


def is_match(query_seq, query_alignment, hit_alignment):
    if query_seq == hit_alignment:
        return True

    query_with_stripped_gaps = query_alignment.strip('-')

    if '-' in query_with_stripped_gaps:
        return False

    hit_with_stripped_gaps = hit_alignment.strip('-')

    if '-' in hit_with_stripped_gaps:
        return False

    return compare_query_and_hit_seqs(query_seq, hit_with_stripped_gaps)


def get_blast_data(pdb_id, chain_id, seq):
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

                    pdb_id = hit_def_parts[0]

                    if pdb_id == '2W9D':
                        print(pdb_id)

                    chain_ids = [x for x in hit_def_parts[2].split(',')]

                    for hsp in hit.find('Hit_hsps'):
                        hsp_qseq = hsp.find('Hsp_qseq').text
                        hsp_hseq = hsp.find('Hsp_hseq').text

                        if not is_match(seq, hsp_qseq, hsp_hseq):
                            continue

                        res.append(BLASTData(pdb_id, chain_ids[0]))

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


def check_names(names):
    if len(list(frozenset(names))) == 1:
        return True

    split_names = list(map(lambda x: x.split(), names))

    if len(list(frozenset(map(lambda x: len(x), split_names)))) != 1:
        return False

    common_set = set([])

    for x in split_names:
        for y in x:
            common_set.add(y.upper())

    return abs(len(list(common_set)) - len(split_names[0])) <= 1


def check_unbound(pdb_id, chain_seqs):
    all_seqs_in_pdb = list(map(lambda x: x[1], fetch_all_sequences(pdb_id)))

    seqs_counts = []

    for chain_seq in chain_seqs:
        seqs_counts.append(0)
        for seq in all_seqs_in_pdb:
            if compare_query_and_hit_seqs(chain_seq, seq):
                seqs_counts[-1] += 1

    # print(pdb_id)
    # print(chain_seqs)
    # print(all_seqs_in_pdb)
    # print(list(map(lambda x: x > 0, seqs_counts)))
    # print(len(retrieve_uniprot_ids(pdb_id)) == 1)

    # we check that for every queried chain there is a matching chain in the
    # given pdb and also we check that given pdb contains only one UniProt
    # structure, what means that it contains only one structure, what means
    # that structure is not in a complex, what means that it's unbound.
    # also if names of all structures in pdb are different in no more than
    # one word (for example, 'my ab heavy chain' and 'my ab light chain)
    # it usually means that structures form one macromolecule,
    # hence their complex is unbound
    return all(map(lambda x: x > 0, seqs_counts)) \
           and (len(retrieve_uniprot_ids(pdb_id)) == 1
                or check_names(retrieve_names(pdb_id)))


def find_unbound_structure(pdb_id, chain_ids, seqs):
    # TODO: идея — искать похожести по количеству совпадающих слов в названии макромолекулы

    # TODO: add memoization to find complexes more effectively
    candidates = [get_blast_data(pdb_id, chain_id, seq) for chain_id, seq in
                  zip(chain_ids, seqs)]

    pdb_ids_in_intersection_prep = reduce(operator.and_,
                                          [set([x.pdb_id for x in candidate])
                                           for
                                           candidate in candidates])

    print(pdb_ids_in_intersection_prep)

    return list(
        filter(lambda x: x.upper() != pdb_id.upper()
                         and check_unbound(x, seqs),
               list(pdb_ids_in_intersection_prep)[:20]))


def find_unbound_conformations(complex):
    unbound_antigen_valid_candidates = \
        find_unbound_structure(complex.pdb_id, complex.antigen_chains,
                               complex.antigen_seqs)

    print(unbound_antigen_valid_candidates)

    unbound_antibody_valid_candidates = \
        find_unbound_structure(complex.pdb_id,
                               [complex.antibody_h_chain,
                                complex.antibody_l_chain],
                               [complex.antibody_h_seq,
                                complex.antibody_l_seq])

    return unbound_antigen_valid_candidates, unbound_antibody_valid_candidates


structures_summary = read_csv('data/sabdab_summary_all.tsv',
                              sep='\t')

# all_complexes = get_bound_complexes(structures_summary)
# load_bound_complexes(all_complexes)

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

comps = get_bound_complexes(structures_summary,
                            list(map(lambda x: x[0], test_structures)))
load_bound_complexes(comps)

for pdb_id, unbound_antibody_id, unbound_antigen_id in test_structures:
    print('processing', pdb_id)

    comp = list(filter(lambda x: x.pdb_id.upper() == pdb_id, comps))[0]
    comp.load_structure()

    unbound_antigen_candidates, unbound_antibody_candidates = \
        find_unbound_conformations(comp)

    print('antigen', 'expected:', unbound_antigen_id, 'got:',
          unbound_antigen_candidates)
    print('antibody', 'expected:', unbound_antibody_id, 'got:',
          unbound_antibody_candidates)

    if unbound_antigen_id not in unbound_antigen_candidates:
        print('MISMATCH! in antigen')

    if unbound_antibody_id not in unbound_antibody_candidates:
        print('MISMATCH! in antibody')

# Проблемы

# 2W9E
# antigen expected: 1QM1 got: ['1HJM', '1QM0', '2LSB', '4N9O', '6DU9', '5YJ5', '1QLX', '1I4M', '1QM3', '1HJN', '1QM2', '2IV5', '1QLZ', '1QM1', '4DGI']
# antibody expected: 2W9D got: []
# MISMATCH! in antibody
#
# 4DN4
# antigen expected: 1DOL got: []
# antibody expected: 4DN3 got: []
# MISMATCH! in antigen
# MISMATCH! in antibody
#
# 4G6M
# antigen expected: 4I1B got: ['9ILB', '7I1B', '2KH2', '1TWM', '5BVP', '2NVH', '5I1B', '6I1B', '1TOO', '1IOB', '4I1B', '2I1B', '1I1B', '5MVZ', '4G6J']
# antibody expected: 4G6K got: []
# MISMATCH! in antibody
#
# 4GXU
# antigen expected: 4I1B got: []
# antibody expected: 4GXV got: []
# MISMATCH! in antigen
# MISMATCH! in antibody

# Запуск 1

# processing 4DN4
# {'4DN4'}
# []
# 2W9D
# 2W9D
# {'4DN3', '2XTJ', '4DN4'}
# antigen expected: 1DOL got: []
# antibody expected: 4DN3 got: ['4DN3']
# MISMATCH! in antigen
# processing 4FQI
# {'4FQI', '6PCX', '6PD6', '3GBM', '6PD5', '4MHH', '3ZP0', '6CFG', '6PD3', '6E3H', '6CF5', '2FK0', '6B3M', '3ZP1'}
# ['6PCX', '6PD6', '3GBM', '6PD5', '3ZP0', '6CFG', '6PD3', '6CF5', '2FK0', '3ZP1']
# 2W9D
# 2W9D
# {'6CNV', '5CJS', '4FQI', '4FQV', '4FQH', '4LLD', '4FQY', '5CJQ'}
# antigen expected: 2FK0 got: ['6PCX', '6PD6', '3GBM', '6PD5', '3ZP0', '6CFG', '6PD3', '6CF5', '2FK0', '3ZP1']
# antibody expected: 4FQH got: ['4FQV', '4FQH', '4FQY']
# processing 4G6J
# {'1I1B', '2KH2', '1ITB', '2I1B', '5MVZ', '4I1B', '7I1B', '3O4O', '9ILB', '4G6J', '4G6M', '1IOB', '5BVP', '5I1B', '2NVH', '4DEP', '6I1B'}
# ['1I1B', '2KH2', '2I1B', '4I1B', '7I1B', '9ILB', '4G6M', '1IOB', '5BVP', '5I1B', '2NVH', '6I1B']
# 2W9D
# 2W9D
# {'4G5Z', '2XTJ', '4G6J'}
# antigen expected: H5N1 got: ['1I1B', '2KH2', '2I1B', '4I1B', '7I1B', '9ILB', '4G6M', '1IOB', '5BVP', '5I1B', '2NVH', '6I1B']
# antibody expected: 4G5Z got: ['4G5Z']
# MISMATCH! in antigen
# processing 4G6M
# {'1I1B', '2KH2', '1ITB', '2I1B', '1TOO', '5MVZ', '4I1B', '7I1B', '4G6J', '9ILB', '3O4O', '4G6M', '1TWM', '1IOB', '5BVP', '5I1B', '2NVH', '4DEP', '6I1B'}
# ['1I1B', '2KH2', '2I1B', '1TOO', '5MVZ', '4I1B', '7I1B', '4G6J', '9ILB', '1TWM', '1IOB', '5BVP', '5I1B', '2NVH', '6I1B']
# 2W9D
# 2W9D
# {'4G6K', '4G6M', '2XTJ'}
# antigen expected: 4I1B got: ['1I1B', '2KH2', '2I1B', '1TOO', '5MVZ', '4I1B', '7I1B', '4G6J', '9ILB', '1TWM', '1IOB', '5BVP', '5I1B', '2NVH', '6I1B']
# antibody expected: 4G6K got: ['4G6K']
# processing 4GXU
# {'3LZF', '1RUZ', '4PY8', '4EEF', '4GXU', '3R2X', '5C0R', '3GBN', '2WRG', '5C0S'}
# ['3LZF', '1RUZ', '4PY8', '4EEF', '3R2X', '3GBN', '2WRG']
# 2W9D
# 2W9D
# {'4LLD', '4GXU', '4GXV'}
# antigen expected: 4I1B got: ['3LZF', '1RUZ', '4PY8', '4EEF', '3R2X', '3GBN', '2WRG']
# antibody expected: 4GXV got: ['4GXV']
# MISMATCH! in antigen

# Запуск 2

# processing 4DN4
# {'4DN4'}
# []
# 2W9D
# 2W9D
# {'2XTJ', '4DN4', '4DN3'}
# antigen expected: 1DOL got: []
# antibody expected: 4DN3 got: ['4DN3']
# MISMATCH! in antigen
# processing 4FQI
# {'6PCX', '6PD5', '2FK0', '6PD6', '3ZP1', '6CF5', '6PD3', '4MHH', '6B3M', '6E3H', '4FQI', '6CFG', '3GBM', '3ZP0'}
# ['6PCX', '6PD5', '2FK0', '6PD6', '3ZP1', '6CF5', '6PD3', '6CFG', '3GBM', '3ZP0']
# 2W9D
# 2W9D
# {'4FQY', '4FQH', '6CNV', '5CJS', '4FQV', '4FQI', '4LLD', '5CJQ'}
# antigen expected: 2FK0 got: ['6PCX', '6PD5', '2FK0', '6PD6', '3ZP1', '6CF5', '6PD3', '6CFG', '3GBM', '3ZP0']
# antibody expected: 4FQH got: ['4FQY', '4FQH', '4FQV']
# processing 4G6J
# {'5I1B', '5BVP', '2KH2', '2I1B', '1IOB', '7I1B', '6I1B', '5MVZ', '4G6M', '4G6J', '1ITB', '4I1B', '4DEP', '2NVH', '1I1B', '3O4O', '9ILB'}
# ['5I1B', '5BVP', '2KH2', '2I1B', '1IOB', '7I1B', '6I1B', '4G6M', '4I1B', '2NVH', '1I1B', '9ILB']
# 2W9D
# 2W9D
# {'2XTJ', '4G6J', '4G5Z'}
# antigen expected: H5N1 got: ['5I1B', '5BVP', '2KH2', '2I1B', '1IOB', '7I1B', '6I1B', '4G6M', '4I1B', '2NVH', '1I1B', '9ILB']
# antibody expected: 4G5Z got: ['4G5Z']
# MISMATCH! in antigen
# processing 4G6M
# {'1TWM', '5I1B', '5BVP', '1TOO', '2KH2', '2I1B', '1IOB', '7I1B', '6I1B', '5MVZ', '4G6M', '4G6J', '1ITB', '4I1B', '4DEP', '2NVH', '1I1B', '3O4O', '9ILB'}
# ['1TWM', '5I1B', '5BVP', '1TOO', '2KH2', '2I1B', '1IOB', '7I1B', '6I1B', '5MVZ', '4G6J', '4I1B', '2NVH', '1I1B', '9ILB']
# 2W9D
# 2W9D
# {'2XTJ', '4G6K', '4G6M'}
# antigen expected: 4I1B got: ['1TWM', '5I1B', '5BVP', '1TOO', '2KH2', '2I1B', '1IOB', '7I1B', '6I1B', '5MVZ', '4G6J', '4I1B', '2NVH', '1I1B', '9ILB']
# antibody expected: 4G6K got: ['4G6K']
# processing 4GXU
# {'3GBN', '3LZF', '1RUZ', '5C0R', '2WRG', '5C0S', '4GXU', '3R2X', '4EEF', '4PY8'}
# ['3GBN', '3LZF', '1RUZ', '2WRG', '3R2X', '4EEF', '4PY8']
# 2W9D
# 2W9D
# {'4GXU', '4LLD', '4GXV'}
# antigen expected: 4I1B got: ['3GBN', '3LZF', '1RUZ', '2WRG', '3R2X', '4EEF', '4PY8']
# antibody expected: 4GXV got: ['4GXV']
# MISMATCH! in antigen
