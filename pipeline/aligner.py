import math
import os
import re
import subprocess

from Bio.PDB import PDBParser, PDBIO
import numpy as np
import pandas as pd

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


def fetch_seqres_and_connect(path):
    with open(path, 'r') as f:
        seqres_lines = []
        connect_lines = []

        for line in f.readlines():
            if line.startswith('SEQRES'):
                seqres_lines.append(line)
            elif line.startswith('CONNECT'):
                connect_lines.append(line)

    return seqres_lines, connect_lines


def get_atoms(struct, only_ca=False):
    atoms = []

    for model in struct:
        for chain in model:
            for residue in chain:
                if only_ca and 'CA' in residue:
                    atoms.append(residue['CA'])
                else:
                    for atom in residue:
                        atoms.append(atom)

    return atoms


def substitute_coords(path, struct):
    atoms_to_coords = {}

    for atom in get_atoms(struct):
        atoms_to_coords[atom.get_serial_number()] = atom.get_coord()

    lines = []

    with open(path, 'r') as f:
        for line in f.readlines():
            if line.startswith('ATOM'):
                serial_number = int(line[6:11])

                x, y, z = atoms_to_coords[serial_number]

                new_coords = '%8.3f%8.3f%8.3f' % (x, y, z)

                line = line[:30] + new_coords + line[30 + len(new_coords):]

            lines.append(line)

    with open(path, 'w') as f:
        f.writelines(lines)


class Candidate:
    def __init__(self, pdb_id_b, ab_pdb_id_u, ag_pdb_id_u, path,
                 ab_chain_ids_b, ag_chain_ids_b,
                 ab_chain_ids_u, ag_chain_ids_u):
        self.pdb_id_b = pdb_id_b
        self.ab_pdb_id_u = ab_pdb_id_u
        self.ag_pdb_id_u = ag_pdb_id_u
        self.path = path
        self.ab_chain_ids_b = ab_chain_ids_b
        self.ag_chain_ids_b = ag_chain_ids_b
        self.ab_chain_ids_u = ab_chain_ids_u
        self.ag_chain_ids_u = ag_chain_ids_u

        self.ab_struct_b = self.fetch_struct('ab_b')
        self.ag_struct_b = self.fetch_struct('ag_b')

        self.ab_path_u = os.path.join(self.path,
                                      self.pdb_id_b + '_ab_u' + DOT_PDB)
        self.ag_path_u = os.path.join(self.path,
                                      self.pdb_id_b + '_ag_u' + DOT_PDB)

        self.ab_struct_u = self.fetch_struct('ab_u')
        self.ag_struct_u = self.fetch_struct('ag_u')

        self.ab_chains_b = Conformation.extract_chains(self.ab_struct_b,
                                                       ab_chain_ids_b)
        self.ag_chains_b = Conformation.extract_chains(self.ag_struct_b,
                                                       ag_chain_ids_b)
        self.ab_chains_u = Conformation.extract_chains(self.ab_struct_u,
                                                       ab_chain_ids_u)
        self.ag_chains_u = Conformation.extract_chains(self.ag_struct_u,
                                                       ag_chain_ids_u)

        self.ab_interface_cas_b, self.ag_interface_cas_b = \
            Conformation.get_interface_atoms(self.ab_struct_b, self.ag_struct_b, only_cas=True)

        self.ab_interface_atoms_u, self.ag_interface_atoms_u = \
            Conformation.get_interface_atoms(self.ab_struct_u, self.ag_struct_u)

    def fetch_struct(self, suff):
        return pdb_parser.get_structure(self.pdb_id_b, os.path.join(self.path,
                                                                    self.pdb_id_b
                                                                    + '_' + suff
                                                                    + DOT_PDB))


def rmsd_by(atoms1, atoms2, by=None):
    atoms1_by = list(filter(lambda x: x.get_id() in by, atoms1)) if by else atoms1
    atoms2_by = list(filter(lambda x: x.get_id() in by, atoms2)) if by else atoms2

    res = 0

    for at1, at2 in zip(atoms1_by, atoms2_by):
        res += np.linalg.norm(at1.coord - at2.coord) ** 2

    res /= len(atoms1_by)

    return round(math.sqrt(res), 3)


def get_rmsds_stats(main_candidate, alternative_candidate):
    main_common_atoms_ab = []
    main_common_atoms_ag = []

    alt_common_atoms_ab = []
    alt_common_atoms_ag = []

    for main_chain, main_chain_id, alt_chain, alt_chain_id in zip(
            main_candidate.ab_chains_u,
            main_candidate.ab_chain_ids_u,
            alternative_candidate.ab_chains_u,
            alternative_candidate.ab_chain_ids_u):
        res = Conformation._matching_atoms_for_chains(main_chain,
                                                      main_candidate.ab_pdb_id_u,
                                                      main_chain_id,
                                                      alt_chain,
                                                      alternative_candidate.ab_pdb_id_u,
                                                      alt_chain_id, only_cas=False)
        main_common_atoms_ab += res[0]
        alt_common_atoms_ab += res[1]

    for main_chain, main_chain_id, alt_chain, alt_chain_id in zip(
            main_candidate.ag_chains_u,
            main_candidate.ag_chain_ids_u,
            alternative_candidate.ag_chains_u,
            alternative_candidate.ag_chain_ids_u):
        res = Conformation._matching_atoms_for_chains(main_chain,
                                                      main_candidate.ag_pdb_id_u,
                                                      main_chain_id,
                                                      alt_chain,
                                                      alternative_candidate.ag_pdb_id_u,
                                                      alt_chain_id, only_cas=False)
        main_common_atoms_ag += res[0]
        alt_common_atoms_ag += res[1]

    main_atoms_tmp = main_common_atoms_ab + main_common_atoms_ag
    alt_atoms_tmp = alt_common_atoms_ab + alt_common_atoms_ag

    main_atoms = []
    alt_atoms = []

    for m_at, a_at in zip(main_atoms_tmp, alt_atoms_tmp):
        if (m_at in main_candidate.ab_interface_atoms_u or
                m_at in main_candidate.ag_interface_atoms_u) and \
            (a_at in alternative_candidate.ab_interface_atoms_u or
                a_at in alternative_candidate.ag_interface_atoms_u):
            main_atoms.append(m_at)
            alt_atoms.append(a_at)

    return rmsd_by(main_atoms, alt_atoms, ['CA']), rmsd_by(main_atoms,
                                                           alt_atoms,
                                                           ['N', 'CA', 'C']), \
           rmsd_by(main_atoms, alt_atoms)


def process_alternative_candidates(path_to_comp, comp_row):
    comp_name = os.path.basename(path_to_comp)

    pdb_id_b, ab_chain_ids_b, ag_chain_ids_b = \
        comp_name_to_pdb_and_chains(comp_name)

    path_to_alternative_candidates = os.path.join(path_to_comp,
                                                  'alternative_candidates')

    if not os.path.exists(path_to_alternative_candidates):
        return

    path_to_alternative_candidates_csv = path_to_alternative_candidates + '.csv'

    alternative_candidates_tbl = pd.read_csv(path_to_alternative_candidates_csv)

    ca_rmsds = []
    ncac_rmsds = []
    all_atoms_rmsds = []

    for candidate_name in alternative_candidates_tbl['candidate_name']:
        try:
            candidate_row = alternative_candidates_tbl[
                alternative_candidates_tbl[
                    'candidate_name'] == candidate_name].iloc[0]

            if candidate_row['type'] != 'U:U':
                continue

            alt_ab_pdb_id_u = candidate_row['ab_pdb_id_u']
            alt_ag_pdb_id_u = candidate_row['ag_pdb_id_u']
            alt_ab_chain_ids_u = candidate_row['ab_chain_ids_u'].split(':')
            alt_ag_chain_ids_u = candidate_row['ag_chain_ids_u'].split(':')

            alt_path_to_comp = os.path.join(path_to_alternative_candidates,
                                            candidate_name)

            for path in os.listdir(alt_path_to_comp):
                pre_path = os.path.join(path_to_comp, path)

                comp_pdb_path = os.path.join(pre_path,
                                             pdb_id_b + DOT_PDB)

                if not os.path.exists(comp_pdb_path):
                    continue

                main_candidate = Candidate(pdb_id_b, comp_row['ab_pdb_id_u'],
                                           comp_row['ag_pdb_id_u'],
                                           pre_path,
                                           ab_chain_ids_b, ag_chain_ids_b,
                                           comp_row['ab_chain_ids_u'].split(
                                               ':'),
                                           comp_row['ag_chain_ids_u'].split(
                                               ':'))

                alt_pre_path = os.path.join(alt_path_to_comp, path)
                alt_pdb_id, _, _ = comp_name_to_pdb_and_chains(candidate_name)
                alternative_candidate = Candidate(alt_pdb_id, alt_ab_pdb_id_u,
                                                  alt_ag_pdb_id_u,
                                                  alt_pre_path,
                                                  alt_ab_pdb_id_u,
                                                  alt_ag_pdb_id_u,
                                                  alt_ab_chain_ids_u,
                                                  alt_ag_chain_ids_u)

                Conformation._inner_align(main_candidate.ab_chain_ids_b,
                                          main_candidate.ab_chains_b,
                                          main_candidate.pdb_id_b,
                                          alternative_candidate.ab_struct_u,
                                          alternative_candidate.ab_chain_ids_u,
                                          alternative_candidate.ab_pdb_id_u,
                                          main_candidate.ab_interface_cas_b)

                substitute_coords(alternative_candidate.ab_path_u,
                                  alternative_candidate.ab_struct_u)

                Conformation._inner_align(main_candidate.ag_chain_ids_b,
                                          main_candidate.ag_chains_b,
                                          main_candidate.pdb_id_b,
                                          alternative_candidate.ag_struct_u,
                                          alternative_candidate.ag_chain_ids_u,
                                          alternative_candidate.ag_pdb_id_u,
                                          main_candidate.ag_interface_cas_b)

                substitute_coords(alternative_candidate.ag_path_u,
                                  alternative_candidate.ag_struct_u)

                if path == 'prepared_schrod':
                    ca_rmsd, ncac_rmsd, all_atoms_rmsd = get_rmsds_stats(main_candidate,
                                                  alternative_candidate)

                    ca_rmsds.append(ca_rmsd)
                    ncac_rmsds.append(ncac_rmsd)
                    all_atoms_rmsds.append(all_atoms_rmsd)
        except Exception as e:
            ca_rmsds.append(None)
            ncac_rmsds.append(None)
            all_atoms_rmsds.append(None)

            print('Can\'t:', candidate_name, e)

    alternative_candidates_tbl['ca_rmsds'] = ca_rmsds
    alternative_candidates_tbl['ncac_rmsds'] = ncac_rmsds
    alternative_candidates_tbl['all_atoms_rmsds'] = all_atoms_rmsds

    os.remove(path_to_alternative_candidates_csv)
    alternative_candidates_tbl.to_csv(path_to_alternative_candidates_csv,
                                      na_rep='NA', index=False)


def process(dir_path):
    counter = 0

    abase_summary_tbl = pd.read_csv(os.path.join(dir_path,
                                                 'abase_summary.csv'))

    data_path = os.path.join(dir_path, 'data')

    for comp_name in os.listdir(data_path):
        pre_path = os.path.join(data_path, comp_name)

        if not os.path.isdir(pre_path):
            continue

        comp_row = abase_summary_tbl[
            abase_summary_tbl[
                'comp_name'] == comp_name].iloc[0]

        counter += 1

        print('Processing {}'.format(comp_name), flush=True)

        process_alternative_candidates(pre_path, comp_row)


if __name__ == '__main__':
    process(os.path.abspath('..'))
