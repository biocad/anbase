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


def get_while_true(curl):
    not_finished = True

    res = None

    while not_finished:
        try:
            res = requests.get(curl)

            if res.content.decode('utf-8'):
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

            if res.content.decode('utf-8'):
                not_finished = False
        except Exception:
            pass

    return res.content.decode('utf-8')


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

    def load_structure(self):
        self.load_structure_from(os.path.join(self.complex_dir_path,
                                              self.db_name + DOT_PDB))

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
            if not seqs:
                print('bad line:', line, 'in', r)
                return fetch_all_sequences(pdb_id)

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

            if is_obsolete(row[PDB_ID]):
                continue

            antigen_chains = row[ANTIGEN_CHAIN].split(' | ')
            complexes.append(Complex(
                row[PDB_ID], sub_nan(row[H_CHAIN]), sub_nan(row[L_CHAIN]),
                antigen_chains,
                sub_nan(row[ANTIGEN_HET_NAME])))

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
    ent_paths = set([])

    with open('could_not_fetch.log', 'w') as could_not_fetch_log:
        pdb_list = PDBList()

        io = PDBIO()

        for comp in complexes:
            pdb_path = os.path.join(comp.complex_dir_path,
                                    comp.db_name + DOT_PDB)

            if os.path.exists(pdb_path):
                if load_structures:
                    comp.load_structure_from(pdb_path)
                print(comp.pdb_id, 'loaded')
                continue

            # if os.path.exists(comp.complex_dir_path):
            #     shutil.rmtree(comp.complex_dir_path)

            if not os.path.exists(comp.complex_dir_path):
                os.mkdir(comp.complex_dir_path)

            ent_path = pdb_list.retrieve_pdb_file(comp.pdb_id,
                                                  file_format='pdb',
                                                  pdir=DB_PATH)

            if not os.path.exists(ent_path):
                print('Not written:', comp.pdb_id)
                print(comp.pdb_id, flush=True, file=could_not_fetch_log)
                continue

            comp.load_structure_from(ent_path)

            needed_chain_ids = [x for x in [comp.antibody_h_chain,
                                            comp.antibody_l_chain] +
                                comp.antigen_chains if x]

            for model in comp.structure:
                for chain in model:
                    if chain.get_id() not in needed_chain_ids:
                        model.detach_child(chain.get_id())

            io.set_structure(comp.structure)
            io.save(pdb_path)

            ent_paths.add(ent_path)

            print(comp.pdb_id, 'loaded')

    for ent_path in ent_paths:
        os.remove(ent_path)


def align_and_check(query_seq, target_seq, pdb_ids, chain_ids, write_log,
                    len_diff):
    cut_off = int(0.05 * len(target_seq))

    alignment = pairwise2.align.localxs(query_seq, target_seq, -10, -10,
                                        penalize_end_gaps=False,
                                        one_alignment_only=True)[0]

    mismatches_count = 0

    query_alignment = alignment[0]
    target_alignment = alignment[1]

    for i in range(len(query_alignment)):
        if query_alignment[i] != '-' and target_alignment[i] != '-' \
                and query_alignment[i] != target_alignment[i]:
            mismatches_count += 1

    if write_log and 0 < mismatches_count <= cut_off:
        with open(MISMATCHED_LOG, 'a') as f:
            f.write(','.join(
                [pdb_ids[0], pdb_ids[1], chain_ids[0].upper(),
                 chain_ids[1].upper(),
                 str(mismatches_count), str(len_diff)]) + '\n')

    return mismatches_count <= cut_off


def compare_query_and_hit_seqs(query_seq, hit_seq, pdb_ids, chain_ids,
                               write_log=False):
    cut_off_half = int(0.1 * len(query_seq) / 2)
    len_diff = abs(len(query_seq) - len(hit_seq))

    if len_diff > 2 * cut_off_half:
        return False

    c1 = align_and_check(query_seq, hit_seq, pdb_ids, chain_ids, write_log,
                         len_diff)

    if c1:
        return True

    c2 = align_and_check(query_seq[cut_off_half:-cut_off_half], hit_seq,
                         pdb_ids, chain_ids, write_log, len_diff)

    if c2:
        return True

    c3 = align_and_check(hit_seq[cut_off_half:-cut_off_half], query_seq,
                         pdb_ids, chain_ids, write_log, len_diff)

    if c3:
        return True

    c4 = align_and_check(query_seq[:-2 * cut_off_half], hit_seq, pdb_ids,
                         chain_ids, write_log, len_diff)

    if c4:
        return True

    c5 = align_and_check(hit_seq[:-2 * cut_off_half], query_seq, pdb_ids,
                         chain_ids, write_log, len_diff)

    if c5:
        return True

    c6 = align_and_check(query_seq[2 * cut_off_half:], hit_seq, pdb_ids,
                         chain_ids, write_log, len_diff)

    if c6:
        return True

    c7 = align_and_check(hit_seq[2 * cut_off_half:], query_seq, pdb_ids,
                         chain_ids, write_log, len_diff)

    return c7


def is_match(query_seq, query_alignment, hit_alignment, pdb_ids, chain_ids):
    if query_seq == hit_alignment:
        return True

    query_with_stripped_gaps = query_alignment.strip('-')

    if '-' in query_with_stripped_gaps:
        return False

    hit_with_stripped_gaps = hit_alignment.strip('-')

    if '-' in hit_with_stripped_gaps:
        return False

    return compare_query_and_hit_seqs(query_seq, hit_with_stripped_gaps,
                                      pdb_ids, chain_ids)


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

                    hit_pdb_id = hit_def_parts[0]

                    hit_chain_ids = [x for x in hit_def_parts[2].split(',')]

                    for hsp in hit.find('Hit_hsps'):
                        hsp_qseq = hsp.find('Hsp_qseq').text
                        hsp_hseq = hsp.find('Hsp_hseq').text

                        if not is_match(seq, hsp_qseq, hsp_hseq,
                                        (pdb_id, hit_pdb_id),
                                        (chain_id, hit_chain_ids[0])):
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
        res.append(pdb.attrib['resolution'])

    return float(res[0])


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


def check_unbound(pdb_id, chain_ids_and_seqs, query_pdb_id):
    all_seqs_in_pdb = fetch_all_sequences(pdb_id)

    chain_matches = defaultdict(list)

    for chain_id, chain_seq in chain_ids_and_seqs:
        for target_chain_id, seq in all_seqs_in_pdb:
            if compare_query_and_hit_seqs(chain_seq, seq,
                                          (query_pdb_id, pdb_id),
                                          (chain_id, target_chain_id),
                                          write_log=True):
                chain_matches[chain_id].append(target_chain_id)

    # we check that for every queried chain there is a matching chain in the
    # given pdb and also we check that given pdb contains only one UniProt
    # structure, what means that it contains only one structure, what means
    # that structure is not in a complex, what means that it's unbound.
    # also if names of all structures in pdb are different in no more than
    # one word (for example, 'my ab heavy chain' and 'my ab light chain)
    # it usually means that structures form one macromolecule,
    # hence their complex is unbound
    c1 = all(map(lambda x: len(chain_matches[x]) > 0, chain_matches.keys()))
    c2 = len(retrieve_uniprot_ids(pdb_id)) == 1
    c3 = check_names(retrieve_names(pdb_id))

    if c1 and (c2 or c3):
        res = []
        for i in range(len(chain_matches[chain_ids_and_seqs[0][0]])):
            res.append(Candidate(pdb_id, list(
                map(lambda x: chain_matches[x][i], chain_matches.keys()))))
        return res

    return []


def sort_and_take_unbound(unbound_candidates):
    unbound_candidates.sort(key=lambda x: -int(x[0]))
    return unbound_candidates[:25]


def find_unbound_structure(pdb_id, chain_ids, seqs):
    candidates = [get_blast_data(pdb_id, chain_id, seq) for chain_id, seq in
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

        res += check_unbound(candidate_id, list(zip(chain_ids, seqs)), pdb_id)

    return res


def sort_and_take_ress(unbound_ress):
    unbound_ress.sort(key=lambda x: retrieve_resolution(x.pdb_id))
    return unbound_ress[:5]


def find_unbound_conformations(complex):
    unbound_antigen_valid_candidates = \
        find_unbound_structure(complex.pdb_id, complex.antigen_chains,
                               complex.antigen_seqs)

    print('unbound antigen:', unbound_antigen_valid_candidates)

    unbound_antibody_valid_candidates = \
        find_unbound_structure(complex.pdb_id,
                               [complex.antibody_h_chain,
                                complex.antibody_l_chain],
                               [complex.antibody_h_seq,
                                complex.antibody_l_seq])

    print('unbound antibody:', unbound_antibody_valid_candidates)

    return sort_and_take_ress(unbound_antigen_valid_candidates), \
           sort_and_take_ress(unbound_antibody_valid_candidates)


structures_summary = read_csv('data/sabdab_summary_all.tsv',
                              sep='\t')


def run_zlab_test():
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

    with open(MISMATCHED_LOG, 'w') as f:
        f.write(
            'bound_id,unbound_id,bound_chain,unbound_chain,mismatches_count,' +
            'len_diff\n')

    comps = get_bound_complexes(structures_summary,
                                list(map(lambda x: x[0], test_structures)))
    load_bound_complexes(comps)

    for pdb_id, unbound_antibody_id, unbound_antigen_id in test_structures:
        print('processing', pdb_id)

        comps_found = list(filter(lambda x: x.pdb_id.upper() == pdb_id, comps))

        for comp in comps_found:
            comp.load_structure()

            unbound_antigen_candidates, unbound_antibody_candidates = \
                find_unbound_conformations(comp)

            print(comp.db_name)

            print('antigen', 'expected:', unbound_antigen_id, 'got:',
                  unbound_antigen_candidates)
            print('antibody', 'expected:', unbound_antibody_id, 'got:',
                  unbound_antibody_candidates)

            if unbound_antigen_id not in list(
                    map(lambda x: x.pdb_id, unbound_antigen_candidates)):
                print('MISMATCH! in antigen')

            if unbound_antibody_id not in list(
                    map(lambda x: x.pdb_id, unbound_antibody_candidates)):
                print('MISMATCH! in antibody')


def collect_unbound_structures():
    comps = get_bound_complexes(structures_summary)
    load_bound_complexes(comps)

    ent_paths = set()

    pdb_list = PDBList()
    pdb_parser = PDBParser()
    io = PDBIO()

    with open('could_not_fetch_final.log', 'w') as could_not_fetch_log:
        for comp in comps:
            comp.load_structure()

            print(comp.db_name)

            unbound_antigen_candidates, unbound_antibody_candidates = \
                find_unbound_conformations(comp)

            def helper_writer(candidates, suf):
                counter = 0
                for candidate in candidates:
                    print('candidate:', candidate)

                    path_to_candidate_pdb = os.path.join(comp.complex_dir_path,
                                                         comp.pdb_id + '_' +
                                                         suf + '_u_' +
                                                         str(counter) + DOT_PDB)

                    ent_path = pdb_list.retrieve_pdb_file(candidate.pdb_id,
                                                          file_format='pdb',
                                                          pdir=
                                                          DB_PATH)

                    if not os.path.exists(ent_path):
                        print('Not written:', comp.pdb_id)
                        print(comp.pdb_id, flush=True,
                              file=could_not_fetch_log)
                        continue

                    structure = pdb_parser.get_structure(candidate.pdb_id,
                                                         ent_path)

                    for model in comp.structure:
                        for chain in model:
                            if chain.get_id() not in candidate.chain_ids:
                                model.detach_child(chain.get_id())

                    io.set_structure(structure)
                    io.save(path_to_candidate_pdb)

                    ent_paths.add(ent_path)

                    counter += 1

            helper_writer(unbound_antigen_candidates, 'AG')
            helper_writer(unbound_antibody_candidates, 'AB')


if __name__ == '__main__':
    if sys.argv and sys.argv[0] == 'test':
        run_zlab_test()
    else:
        collect_unbound_structures()
