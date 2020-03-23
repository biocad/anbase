import os
import re
import subprocess

from Bio.PDB import PDBParser, PDBIO

from process_unbound_data import Conformation

DOT_PDB = '.pdb'
DOT_MAE = '.mae'

CHAINS_SEPARATOR = '+'

pdb_parser = PDBParser()
pdb_io = PDBIO()


def comp_name_to_pdb_and_chains(comp_name):
    [pdb_id, chains] = comp_name.split('_')[:2]
    ab_chains_s, ag_chains_s = chains.split('|')

    ab_chains = ab_chains_s.split(CHAINS_SEPARATOR)
    ag_chains = ag_chains_s.split(CHAINS_SEPARATOR)

    return pdb_id, ab_chains, ag_chains


def filter_chains(struct, chain_ids):
    for model in struct:
        chains = list(model.get_chains())
        for chain in chains:
            if chain.get_id() not in chain_ids:
                model.detach_child(chain.get_id())


def pdb_to_mae(pdb_path):
    path_to_mae = pdb_path[:-4] + DOT_MAE

    if os.path.exists(path_to_mae):
        return

    command = '$SCHRODINGER/utilities/structconvert ' \
              '-ipdb \'{}\' -omae \'{}\''.format(pdb_path, path_to_mae)
    subprocess.call(command, stdout=subprocess.PIPE, shell=True)


def process(dir_path):
    def accept_path(p):
        found = re.search('(...._(.\+.|.)\|(.\+.\+.\+.|.\+.\+.|.\+.|.))', p)
        return found and found.endpos == len(p)

    l = list(os.walk(dir_path))

    counter = 0
    for root, dirnames, _ in l:
        counter += 1

        print('Processing [{}/{}]:'.format(counter, len(l)), root, flush=True)

        for file in dirnames:
            if accept_path(os.path.basename(root)):
                comp_name = os.path.basename(root)
                pdb_id, ab_chain_ids, ag_chain_ids = \
                    comp_name_to_pdb_and_chains(comp_name)

                pre_path = os.path.join(root, file)

                comp_pdb_path = os.path.join(pre_path, pdb_id + DOT_PDB)

                if not os.path.exists(comp_pdb_path):
                    continue

                struct = pdb_parser.get_structure(pdb_id,
                                                  os.path.join(pre_path,
                                                               pdb_id +
                                                               DOT_PDB))

                complex_ab_b_path = os.path.join(pre_path,
                                                 pdb_id + '_ab_b' + DOT_PDB)
                complex_ag_b_path = os.path.join(pre_path,
                                                 pdb_id + '_ag_b' + DOT_PDB)

                ab_struct_b = struct.copy()
                filter_chains(ab_struct_b, ab_chain_ids)

                Conformation.write_structure(ab_struct_b, complex_ab_b_path,
                                             pdb_id,
                                             {x: x for x in ab_chain_ids})

                if file.endswith('prepared_schrod'):
                    pdb_to_mae(complex_ab_b_path)

                ag_struct_b = struct.copy()
                filter_chains(ag_struct_b, ag_chain_ids)

                Conformation.write_structure(ag_struct_b, complex_ag_b_path,
                                             pdb_id,
                                             {x: x for x in ag_chain_ids})

                if file.endswith('prepared_schrod'):
                    pdb_to_mae(complex_ag_b_path)


if __name__ == '__main__':
    process(os.path.abspath('../data'))
