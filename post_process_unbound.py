import os
from collections import defaultdict
import pandas as pd
from Bio import pairwise2
from Bio.PDB import PDBParser, Superimposer
from Bio.PDB.Polypeptide import dindex_to_1, d3_to_index, PPBuilder
import numpy as np

from collect_db import fetch_all_sequences, AG, AB, DB_PATH, DOT_PDB, \
    retrieve_pdb


class Conformation:
    pdb_parser = PDBParser()
    super_imposer = Superimposer()
    peptides_builder = PPBuilder()

    def __init__(self, pdb_id_b, heavy_chain_id_b,
                 light_chain_id_b, antigen_chain_ids_b,
                 pdb_ab_id_u, heavy_chain_id_u, light_chain_id_u,
                 pdb_ag_id_u, antigen_chain_ids_u, is_ab_u, is_ag_u):
        self.pdb_id_b = pdb_id_b
        self.heavy_chain_id_b = heavy_chain_id_b
        self.light_chain_id_b = light_chain_id_b
        self.antigen_chain_ids_b = antigen_chain_ids_b
        self.pdb_ab_id_u = pdb_ab_id_u
        self.heavy_chain_id_u = heavy_chain_id_u
        self.light_chain_id_u = light_chain_id_u
        self.pdb_ag_id_u = pdb_ag_id_u
        self.antigen_chain_ids_u = antigen_chain_ids_u

        self.complex_structure_b = self._load_structure(pdb_id_b)
        self.ab_structure_u = self._load_structure(pdb_ab_id_u)
        self.ag_structure_u = self._load_structure(pdb_ag_id_u)

        self.is_ab_u = is_ab_u
        self.is_ag_u = is_ag_u

        self.ab_chains_b = self.extract_chains(self.complex_structure_b,
                                               [self.heavy_chain_id_b,
                                                self.light_chain_id_b])
        self.ag_chains_b = self.extract_chains(self.complex_structure_b,
                                               self.antigen_chain_ids_b)

        self.ab_atoms_b = []
        self.ag_atoms_b = []

        for chain in self.ab_chains_b:
            self.ab_atoms_b += self.extract_cas(chain)

        for chain in self.ag_chains_b:
            self.ag_atoms_b += self.extract_cas(chain)

        self.ab_interface_cas, self.ag_interface_cas = self.get_interface_cas()

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
    def _matching_atoms_for_chains(chain1, chain2):
        #TODO: тут
        l1_seq = Conformation.peptides_builder.build_peptides(chain1)[0].get_sequence()
        l2_seq = Conformation.peptides_builder.build_peptides(chain2)[0].get_sequence()

        alignment = pairwise2.align.localxs(l2_seq, l1_seq, -5, -1,
                                            penalize_end_gaps=False,
                                            one_alignment_only=True)[0]

        a = 42

    def align_ab(self):
        if not self.is_ab_u:
            return

        [heavy_chain_u, light_chain_u] = self.extract_chains(
            self.ab_structure_u, [self.heavy_chain_id_u,
                                  self.light_chain_id_u])

        ref_atoms = self.ab_atoms_b
        sample_atoms = self.extract_cas(heavy_chain_u) + \
                       self.extract_cas(light_chain_u)

        self._align_two_lists_of_atoms(ref_atoms, sample_atoms)
        self.super_imposer.apply(self.ab_structure_u.get_atoms())

        print(self.super_imposer.rms)

    def align_ag(self):
        if not self.is_ag_u:
            return

        ag_chains_u = self.extract_chains(self.ag_structure_u,
                                          self.antigen_chain_ids_u)

        ref_atoms = self.ag_atoms_b
        sample_atoms = []

        for chain in ag_chains_u:
            sample_atoms += self.extract_cas(chain)

        self._matching_atoms_for_chains(self.ag_chains_b[0], ag_chains_u[0])
        self.super_imposer.apply(self.ab_structure_u.get_atoms())

        print(self.super_imposer.rms)


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

    for pdb_ag_id_u, chainss_ag in ag_pdbs_with_chains:
        for chains_ag in chainss_ag:
            chains_ag_split = chains_ag.split(':')
            for pdb_ab_id_u, chainss_ab in ab_pdbs_with_chains:
                for chains_ab in chainss_ab:
                    [heavy_chain_id_u, light_chain_id_u] = chains_ab.split(':')
                    conformation = Conformation(pdb_id_b, heavy_chain_id_b,
                                                light_chain_id_b,
                                                ag_chain_ids_b, pdb_ab_id_u,
                                                heavy_chain_id_u,
                                                light_chain_id_u,
                                                pdb_ag_id_u, chains_ag_split,
                                                is_ab_u, is_ag_u)
                    conformation.align_ab()
                    conformation.align_ag()
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
