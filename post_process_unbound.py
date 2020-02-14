import Bio
import os
from collections import defaultdict
import pandas as pd
from Bio import pairwise2
from Bio.PDB import PDBParser, Superimposer, Chain, PDBIO
from Bio.PDB.Polypeptide import dindex_to_1, d3_to_index, PPBuilder
import numpy as np
from Bio.PDB.StructureBuilder import StructureBuilder

from collect_db import fetch_all_sequences, AG, AB, DB_PATH, DOT_PDB, \
    retrieve_pdb, fetch_sequence


def get_unpacked_list(self):
    """
    Returns all atoms from the residue,
    in case of disordered, keep only first alt loc and remove the alt-loc tag
    """
    atom_list = self.get_list()
    undisordered_atom_list = []
    for atom in atom_list:
        if atom.is_disordered():
            atom.altloc = " "
            undisordered_atom_list.append(atom)
        else:
            undisordered_atom_list.append(atom)
    return undisordered_atom_list


Bio.PDB.Residue.Residue.get_unpacked_list = get_unpacked_list


class Conformation:
    pdb_parser = PDBParser()
    super_imposer = Superimposer()
    peptides_builder = PPBuilder()
    pdb_io = PDBIO()

    def __init__(self, pdb_id_b, heavy_chain_id_b,
                 light_chain_id_b, ag_chain_ids_b,
                 ab_pdb_id_u, heavy_chain_id_u, light_chain_id_u,
                 ag_pdb_id_u, ag_chain_ids_u, is_ab_u, is_ag_u, candidate_id):
        self.pdb_id_b = pdb_id_b
        self.heavy_chain_id_b = heavy_chain_id_b
        self.light_chain_id_b = light_chain_id_b
        self.ag_chain_ids_b = ag_chain_ids_b
        self.ab_pdb_id_u = ab_pdb_id_u
        self.heavy_chain_id_u = heavy_chain_id_u
        self.light_chain_id_u = light_chain_id_u
        self.ag_pdb_id_u = ag_pdb_id_u
        self.ag_chain_ids_u = ag_chain_ids_u

        self.complex_structure_b = self._load_structure(pdb_id_b)
        self.ab_structure_u = self._load_structure(ab_pdb_id_u)
        self.ag_structure_u = self._load_structure(ag_pdb_id_u)

        self.is_ab_u = is_ab_u
        self.is_ag_u = is_ag_u

        self.ab_chains_b = self.extract_chains(self.complex_structure_b,
                                               [self.heavy_chain_id_b,
                                                self.light_chain_id_b])
        self.ag_chains_b = self.extract_chains(self.complex_structure_b,
                                               self.ag_chain_ids_b)

        self.ab_atoms_b = []
        self.ag_atoms_b = []

        for chain in self.ab_chains_b:
            self.ab_atoms_b += self.extract_cas(chain)

        for chain in self.ag_chains_b:
            self.ag_atoms_b += self.extract_cas(chain)

        self.ab_interface_cas, self.ag_interface_cas = self.get_interface_cas()

        self.candidate_id = candidate_id

    @staticmethod
    def extract_chains(structure, chain_ids):
        chains = []

        for chain in structure:
            if chain.get_id() in chain_ids:
                chains.append(chain)

        return chains

    @staticmethod
    def extract_cas(chain):
        cas = []

        for res in chain:
            if 'CA' in res:
                cas.append(res['CA'])

        return cas

    def get_interface_cas(self):
        interface_cutoff = 10

        ab_interface_cas = []
        ag_interface_cas = []

        for ab_at in self.ab_atoms_b:
            for ag_at in self.ag_atoms_b:
                if np.linalg.norm(
                        ab_at.coord - ag_at.coord) < interface_cutoff:
                    ab_interface_cas.append(ab_at)
                    ag_interface_cas.append(ag_at)

        return ab_interface_cas, ag_interface_cas

    @staticmethod
    def _load_structure(pdb_id):
        return \
            Conformation.pdb_parser.get_structure(pdb_id,
                                                  retrieve_pdb(pdb_id))[0]

    @staticmethod
    def _matching_atoms_for_chains(chain1, pdb_id1, chain_id1, chain2, pdb_id2,
                                   chain_id2):
        def extract_seq(chain):
            seq = ''

            for x in Conformation.peptides_builder.build_peptides(chain):
                seq += str(x.get_sequence())

            return seq

        def extract_peps(chain):
            peps = []

            for x in Conformation.peptides_builder.build_peptides(chain):
                peps += x

            return peps

        def get_ids_from_chain(chain, seq, ids_in_seq):
            struct_seq = extract_seq(chain)

            alignment_loc = \
                pairwise2.align.localxs(struct_seq, seq, -5, -1,
                                        penalize_end_gaps=False,
                                        one_alignment_only=True)[0]

            counter = -1
            counter_seq = -1

            res = []

            for i in range(len(alignment_loc[0])):
                if alignment_loc[0][i] == '-' and alignment_loc[1][i] == '-':
                    continue
                elif alignment_loc[0][i] == '-':
                    counter_seq += 1
                    continue
                elif alignment_loc[1][i] == '-':
                    counter += 1
                    continue
                else:
                    counter += 1
                    counter_seq += 1

                if counter_seq in ids_in_seq:
                    res.append((counter_seq, counter))

            return {key: value for key, value in res}

        seq1 = fetch_sequence(pdb_id1, chain_id1)
        seq2 = fetch_sequence(pdb_id2, chain_id2)

        alignment = \
            pairwise2.align.localxs(seq1, seq2, -5, -1,
                                    penalize_end_gaps=False,
                                    one_alignment_only=True)[0]

        counter1 = -1
        counter2 = -1

        ids_in_seq1 = []
        ids_in_seq2 = []

        for i in range(len(alignment[0])):
            if alignment[0][i] == '-' and alignment[1][i] == '-':
                continue
            elif alignment[0][i] == '-':
                counter2 += 1
                continue
            elif alignment[1][i] == '-':
                counter1 += 1
                continue
            else:
                counter1 += 1
                counter2 += 1

            ids_in_seq1.append(counter1)
            ids_in_seq2.append(counter2)

        peps1 = extract_peps(chain1)
        peps2 = extract_peps(chain2)

        ids1 = get_ids_from_chain(chain1, seq1, ids_in_seq1)
        ids2 = get_ids_from_chain(chain2, seq2, ids_in_seq2)

        mutual_ids = frozenset(ids1.keys()) & frozenset(ids2.keys())

        atoms1 = [peps1[ids1[i]]['CA'] for i in mutual_ids]
        atoms2 = [peps2[ids2[i]]['CA'] for i in mutual_ids]

        return atoms1, atoms2

    def align_ab(self):
        if not self.is_ab_u:
            return

        [heavy_chain_u, light_chain_u] = self.extract_chains(
            self.ab_structure_u, [self.heavy_chain_id_u,
                                  self.light_chain_id_u])

        heavy_atoms1, heavy_atoms2 = self._matching_atoms_for_chains(
            self.ab_chains_b[0],
            self.pdb_id_b,
            self.heavy_chain_id_b,
            heavy_chain_u,
            self.ab_pdb_id_u,
            self.heavy_chain_id_u)

        light_atoms1, light_atoms2 = self._matching_atoms_for_chains(
            self.ab_chains_b[1],
            self.pdb_id_b,
            self.light_chain_id_b,
            light_chain_u,
            self.ab_pdb_id_u,
            self.light_chain_id_u)

        self.super_imposer.set_atoms(heavy_atoms1 + light_atoms1,
                                     heavy_atoms2 + light_atoms2)
        self.super_imposer.apply(self.ab_structure_u.get_atoms())

        print(self.super_imposer.rms)

    def align_ag(self):
        if not self.is_ag_u:
            return

        ag_chains_u = self.extract_chains(self.ag_structure_u,
                                          self.ag_chain_ids_u)

        atoms1 = []
        atoms2 = []

        for i in range(len(ag_chains_u)):
            tmp_atoms1, tmp_atoms2 = self._matching_atoms_for_chains(
                self.ag_chains_b[i],
                self.pdb_id_b,
                self.ag_chain_ids_b[i],
                ag_chains_u[i],
                self.ag_pdb_id_u,
                self.ag_chain_ids_u[i])

            atoms1 += tmp_atoms1
            atoms2 += tmp_atoms2

        self.super_imposer.set_atoms(atoms1, atoms2)
        self.super_imposer.apply(self.ag_structure_u.get_atoms())

        print(self.super_imposer.rms)

    def write_candidate(self):
        path = os.path.join(self.pdb_id_b, str(self.candidate_id))
        name_prefix = os.path.join(path, self.ab_pdb_id_u + '_' + self.ag_pdb_id_u)

        if not os.path.exists(path):
            os.makedirs(path)

        sb = StructureBuilder()

        sb.init_structure('complex')
        sb.init_model(0)

        for chain in self.ab_structure_u.copy():
            sb.model.add(chain)

        for chain in self.ag_structure_u.copy():
            sb.model.add(chain)

        self.pdb_io.set_structure(sb.structure)
        self.pdb_io.save(
            name_prefix + '_complex' + ('_u' if self.is_ab_u else '_b')
            + DOT_PDB)

        self.pdb_io.set_structure(self.ab_structure_u)
        self.pdb_io.save(
            name_prefix + '_r' + ('_u' if self.is_ab_u else '_b')
            + DOT_PDB)

        self.pdb_io.set_structure(self.ag_structure_u)
        self.pdb_io.save(
            name_prefix + '_l' + ('_u' if self.is_ag_u else '_b')
            + DOT_PDB)


def process_csv(csv):
    data = defaultdict(list)

    for i in range(len(csv)):
        data[csv.iloc[i]['db_name']].append((csv.iloc[i]['type'],
                                             csv.iloc[i]['candidate_pdb_id'],
                                             csv.iloc[i][
                                                 'candidate_chain_names']))

    return data


def get_pbds_with_chains(candidates, ty):
    pdbs_to_chains = defaultdict(list)

    for x in candidates:
        if x[0] != ty:
            continue

        pdbs_to_chains[x[1]].append(x[2])

    return list(pdbs_to_chains.items())


def process_candidates(db_name, candidates):
    db_name_split = db_name.split('_')

    pdb_id_b = db_name_split[0]
    chains_sep = list(map(lambda x: db_name_split[1][x],
                          range(len(db_name_split[1]))))

    # TODO: NO VHHs?
    heavy_chain_id_b = chains_sep[0]
    light_chain_id_b = chains_sep[1]
    ag_chain_ids_b = chains_sep[2:]

    ag_pdbs_with_chains = get_pbds_with_chains(candidates, AG)
    ab_pdbs_with_chains = get_pbds_with_chains(candidates, AB)

    is_ab_u = True
    is_ag_u = True

    if not ab_pdbs_with_chains:
        is_ab_u = False
        ab_pdbs_with_chains = [(pdb_id_b, [heavy_chain_id_b + ':' +
                                           light_chain_id_b])]

    if not ag_pdbs_with_chains:
        is_ag_u = False
        ag_pdbs_with_chains = [(pdb_id_b, ag_chain_ids_b)]

    res = []

    counter = -1

    for ag_pdb_id_u, chainss_ag in ag_pdbs_with_chains:
        for chains_ag in chainss_ag:
            chains_ag_split = chains_ag.split(':')
            for ab_pdb_id_u, chainss_ab in ab_pdbs_with_chains:
                for chains_ab in chainss_ab:
                    counter += 1

                    [heavy_chain_id_u, light_chain_id_u] = chains_ab.split(':')
                    conformation = Conformation(pdb_id_b, heavy_chain_id_b,
                                                light_chain_id_b,
                                                ag_chain_ids_b, ab_pdb_id_u,
                                                heavy_chain_id_u,
                                                light_chain_id_u,
                                                ag_pdb_id_u, chains_ag_split,
                                                is_ab_u, is_ag_u, counter)
                    conformation.align_ab()
                    conformation.align_ag()
                    conformation.write_candidate()
                    res.append(conformation)

    return res

    # 1. Если гомомер, то дотаскиваем цепи
    # 2. Если есть молекула вблизи интефрейса взаимодействия, то помечаем
    # 3. Удаляем HETATOMы


def process_unbound(path_to_unbound_csv):
    prepared = process_csv(pd.read_csv(path_to_unbound_csv))

    for key, value in prepared.items():
        process_candidates(key, value)


process_unbound('unbound_data.csv')
