import math
import requests

from Bio.PDB import PDBList, PDBParser, PDBIO, Selection, Polypeptide
from pandas import read_csv
import os
import shutil
from xml.etree import ElementTree
from Bio.PDB.Polypeptide import PPBuilder
import response
import collections
import Bio.PDB

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


class SAbDabEntry:
    pdb_parser = PDBParser()

    def __init__(self, pdb_id, h_chain, l_chain, antigen_chain,
                 antigen_het_name):
        self.pdb_id = pdb_id
        self.h_chain = h_chain
        self.l_chain = l_chain
        self.antigen_chains = antigen_chain
        self.antigen_het_name = antigen_het_name
        self.structure = None

        self.antigen_seqs = None
        self.antibody_h_seq = None
        self.antibody_l_seq = None

    def load_structure(self, path):
        self.structure = self.pdb_parser.get_structure(self.pdb_id, path)

        chains_list = list(self.structure.get_chains())
        ppb = PPBuilder()

        if self.antigen_chains is not None:
            self.antigen_seqs = []

            for antigen_chain in self.antigen_chains:
                chain = \
                    [x for x in chains_list if x.get_id() == antigen_chain][0]
                self.antigen_seqs.append(
                    ppb.build_peptides(chain)[0].get_sequence())

        if self.h_chain is not None:
            chain = [x for x in chains_list if x.get_id() == self.h_chain][0]
            self.antibody_h_seq = ppb.build_peptides(chain)[0].get_sequence()

        if self.h_chain is not None:
            chain = [x for x in chains_list if x.get_id() == self.l_chain][0]
            self.antibody_l_seq = ppb.build_peptides(chain)[0].get_sequence()


def get_bound_complexes(sabdab_summary_df):
    def sub_nan(val):
        if isinstance(val, float) and math.isnan(val):
            return None
        return val

    complexes = []

    for _, row in sabdab_summary_df.iterrows():
        # TODO: remove 6nyq
        if sub_nan(row[ANTIGEN_TYPE]) is not None and row[PDB_ID] == '6nyq':
            antigen_chains = row[ANTIGEN_CHAIN].split(' | ')
            complexes.append(SAbDabEntry(
                row[PDB_ID], sub_nan(row[H_CHAIN]), sub_nan(row[L_CHAIN]),
                antigen_chains,
                sub_nan(row[ANTIGEN_HET_NAME])))

    return complexes


class BLASTData:
    def __init__(self, pdb_id, chain_id):
        self.pdb_id = pdb_id
        self.chain_id = chain_id


def load_bound_complexes(complexes, load_structures=False):
    with open('could_not_fetch.log', 'w') as could_not_fetch_log:
        pdb_list = PDBList()

        io = PDBIO()

        for comp in complexes:
            pdb_dir_path = os.path.join(DB_PATH, comp.pdb_id)
            pdb_path = os.path.join(pdb_dir_path, comp.pdb_id + DOT_PDB)

            if os.path.exists(pdb_path):
                if load_structures:
                    comp.load_structure(pdb_path)
                print(comp.pdb_id, 'loaded')
                continue

            if os.path.exists(pdb_dir_path):
                shutil.rmtree(pdb_dir_path)

            os.mkdir(pdb_dir_path)

            ent_path = pdb_list.retrieve_pdb_file(comp.pdb_id,
                                                  file_format='pdb',
                                                  pdir=DB_PATH)

            if not os.path.exists(ent_path):
                print('Not written:', comp.pdb_id)
                print(comp.pdb_id, flush=True, file=could_not_fetch_log)
                continue

            comp.load_structure(ent_path)

            needed_chain_ids = [x for x in [comp.h_chain, comp.l_chain] +
                                comp.antigen_chain if x is not None]

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
                    structure_id = int(hit_def_parts[1])
                    chain_ids = [x for x in hit_def_parts[2].split(',')]

                    # TODO: structure_id != 1 — good decision?
                    if len(chain_ids) != 1 or structure_id != 1:
                        continue

                    for hsp in hit:
                        # TODO: error here
                        hsp_qseq = hsp.find('Hsp_qseq').text

                        # TODO: add unbound check?
                        if len(hsp_qseq) != len(seq):
                            continue

                        res.append(BLASTData(pdb_id, chain_ids[1]))

    return res


def find_unbound_conformations(complex):
    return None


structures_summary = read_csv('data/sabdab_summary_all.tsv',
                              sep='\t')

# all_complexes = get_bound_complexes(structures_summary)
# load_bound_complexes(all_complexes)

comp = get_bound_complexes(structures_summary)[0]
comp.load_structure(os.path.join(DB_PATH, '6nyq', '6nyq.pdb'))
print(get_blast_data('6nyq', 'H', comp.antibody_h_seq))
