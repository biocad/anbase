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

        with open(fasta_path, 'w') as f:
            f.write(fasta[0] + '\n' + fasta[1])

        return fasta[1]


def fetch_sequence(pdb_id, chain_id):
    url = 'https://www.rcsb.org/pdb/download/downloadFastaFiles.do'
    r = requests.post(url, {'structureIdList': pdb_id,
                            'compressionType': 'uncompressed'})

    writing = False
    seq = ''

    for x in r.content.decode('utf-8').split():
        if ':' + chain_id + '|' in x:
            writing = True
            continue
        elif writing and '|' in x:
            writing = False

        if writing:
            seq += x

    return seq


def get_bound_complexes(sabdab_summary_df):
    def sub_nan(val):
        if isinstance(val, float) and math.isnan(val):
            return None
        return val

    complexes = []

    for _, row in sabdab_summary_df.iterrows():
        # if antigen's type is in lower case, it means that antigen is no good
        # for us, because it's a small molecule
        if sub_nan(row[ANTIGEN_TYPE]) and row[ANTIGEN_TYPE].islower() and row[
            PDB_ID] == '1dqj':
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

            needed_chain_ids = [x for x in [comp.h_chain, comp.l_chain] +
                                comp.antigen_chain if x]

            for model in comp.structure:
                for chain in model:
                    if chain.get_id() not in needed_chain_ids:
                        model.detach_child(chain.get_id())

            io.set_structure(comp.structure)
            io.save(pdb_path)

            os.remove(ent_path)

            print(comp.pdb_id, 'loaded')


def get_blast_data(pdb_id, chain_id, seq):
    curl = 'https://www.rcsb.org/pdb/rest/getBlastPDB2?structureId' \
           '={}&chainId={}&eCutOff=10.0&matrix=BLOSUM62&outputFormat=XML'. \
        format(pdb_id, chain_id)

    r = requests.get(curl)
    xml = ElementTree.fromstring(r.content)

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
                    print(pdb_id)
                    chain_ids = [x for x in hit_def_parts[2].split(',')]

                    for hsp in hit.find('Hit_hsps'):
                        hsp_qseq = hsp.find('Hsp_qseq').text
                        hsp_hseq = hsp.find('Hsp_hseq').text

                        if pdb_id == '1DQQ':
                            print('seq', seq)
                            print(hsp_qseq)
                            print(hsp_hseq)

                        if hsp_hseq != seq:
                            continue

                        res.append(BLASTData(pdb_id, chain_ids[0]))

    return res


def find_unbound_structure(pdb_id, chain_ids, seqs):
    candidates = [get_blast_data(pdb_id, chain_id, seq) for chain_id, seq in
                  zip(chain_ids, seqs)]

    pdb_ids_in_intersection_prep = reduce(operator.and_,
                                          [set([x.pdb_id for x in candidate])
                                           for
                                           candidate in candidates])

    print(pdb_ids_in_intersection_prep)


def find_unbound_conformations(complex):
    return None


structures_summary = read_csv('data/sabdab_summary_all.tsv',
                              sep='\t')

# all_complexes = get_bound_complexes(structures_summary)
# load_bound_complexes(all_complexes)

comp = get_bound_complexes(structures_summary)[0]
comp.load_structure()
find_unbound_structure(comp.pdb_id,
                           [comp.antibody_h_chain, comp.antibody_l_chain],
                           [comp.antibody_h_seq, comp.antibody_l_seq])

# print(a)
