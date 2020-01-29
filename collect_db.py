import math

from Bio.PDB import PDBList, PDBParser
from pandas import read_csv
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


class SAbDabEntry:
    def __init__(self, pdb_id, h_chain, l_chain, antigen_chain,
                 antigen_het_name):
        self.pdb_id = pdb_id
        self.h_chain = h_chain
        self.l_chain = l_chain
        self.antigen_chain = antigen_chain
        self.antigen_het_name = antigen_het_name
        self.structure = None


def get_bound_complexes(sabdab_summary_df):
    def sub_nan(val):
        if isinstance(val, float) and math.isnan(val):
            return None
        return val

    complexes = []

    for _, row in sabdab_summary_df.iterrows():
        if sub_nan(row[ANTIGEN_TYPE]) is not None:
            antigen_chains = row[ANTIGEN_CHAIN].split(' | ')
            complexes.append(SAbDabEntry(
                row[PDB_ID], sub_nan(row[H_CHAIN]), sub_nan(row[L_CHAIN]),
                antigen_chains,
                sub_nan(row[ANTIGEN_HET_NAME])))

    return complexes


def load_bound_complexes(complexes):
    pdbl = PDBParser()

    for comp in complexes:
        comp.structure = pdbl.get_structure(comp.pdb_id,
                                            'resources/all_structures/raw/'
                                            + comp.pdb_id + '.pdb')[0]


structures_summary = read_csv('resources/sabdab_summary_all.tsv',
                              sep='\t')

complexes = get_bound_complexes(structures_summary)
load_bound_complexes(complexes)
