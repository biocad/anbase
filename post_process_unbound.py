import string
from xml.etree import ElementTree

import Bio
import os
from collections import defaultdict
import pandas as pd
from Bio import pairwise2
from Bio.PDB import PDBParser, Superimposer, Chain, PDBIO, Select
from Bio.PDB.Polypeptide import dindex_to_1, d3_to_index, PPBuilder
import numpy as np
from Bio.PDB.StructureBuilder import StructureBuilder

from collect_db import fetch_all_sequences, AG, AB, DB_PATH, DOT_PDB, \
    retrieve_pdb, fetch_sequence, with_timeout, memoize, get_while_true, \
    compare_query_and_hit_seqs

# copy=Bio.PDB.Atom.copy
# def myCopy(self):
#     shallow = copy.copy(self)
#     for child in self.child_dict.values():
#         shallow.disordered_add(child.copy())
#     return shallow
# Bio.PDB.Atom.DisorderedAtom.copy=myCopy

# def get_unpacked_list(self):
#     """
#     Returns all atoms from the residue,
#     in case of disordered, keep only first alt loc and remove the alt-loc tag
#     """
#     atom_list = self.get_list()
#     undisordered_atom_list = []
#     for atom in atom_list:
#         if atom.is_disordered():
#             atom.altloc = " "
#             undisordered_atom_list.append(atom)
#         else:
#             undisordered_atom_list.append(atom)
#     return undisordered_atom_list
#
#
# Bio.PDB.Residue.Residue.get_unpacked_list = get_unpacked_list

FILTERED_STRUCTURES_CSV = 'filtered_for_unboundness.csv'
REJECTED_STRUCTURES_CSV = 'rejected_for_unboundness.csv'

FILTERED_COMPLEXES_CSV = 'filtered_complexes.csv'
REJECTED_COMPLEXES_CSV = 'rejected_complexes.csv'


class NotDisordered(Select):
    # this crutch is needed due to the fact that biopython is bad at handling
    # atoms with alternate locations. So we just delete them
    def accept_atom(self, atom):
        if not atom.is_disordered() or atom.get_altloc() == 'A':
            if atom.get_altloc() == 'A':
                atom.altloc = ' '
            return True
        return False


class Conformation:
    pdb_parser = PDBParser()
    super_imposer = Superimposer()
    peptides_builder = PPBuilder()
    pdb_io = PDBIO()

    def __init__(self, pdb_id_b, assembly_id_b, heavy_chain_id_b,
                 light_chain_id_b, ag_chain_ids_b,
                 ab_pdb_id_u, ab_assembly_id, heavy_chain_id_u,
                 light_chain_id_u,
                 ag_pdb_id_u, ag_assembly_id, ag_chain_ids_u, is_ab_u, is_ag_u,
                 candidate_id):
        self.pdb_id_b = pdb_id_b
        self.assembly_id_b = assembly_id_b
        self.heavy_chain_id_b = heavy_chain_id_b
        self.light_chain_id_b = light_chain_id_b
        self.ag_chain_ids_b = ag_chain_ids_b
        self.ab_pdb_id_u = ab_pdb_id_u
        self.heavy_chain_id_u = heavy_chain_id_u
        self.light_chain_id_u = light_chain_id_u
        self.ag_pdb_id_u = ag_pdb_id_u
        self.ag_chain_ids_u = ag_chain_ids_u
        self.ab_assembly_id = ab_assembly_id
        self.ag_assembly_id = ag_assembly_id

        self.is_ab_u = is_ab_u

        self.is_ag_u = is_ag_u

        self.complex_structure_b = self._load_structure(pdb_id_b,
                                                        self.assembly_id_b)
        if self.is_ab_u:
            self.ab_structure_u = self._load_structure(ab_pdb_id_u,
                                                       self.ab_assembly_id)
        else:
            self.ab_structure_u = self.complex_structure_b.copy()

            for model in self.ab_structure_u:
                for chain in model:
                    if chain.get_id() not in [self.heavy_chain_id_b,
                                              self.light_chain_id_b]:
                        model.detach_child(chain.get_id())

        if self.is_ag_u:
            self.ag_structure_u = self._load_structure(ag_pdb_id_u,
                                                       self.ag_assembly_id)
        else:
            self.ag_structure_u = self.complex_structure_b.copy()

            for model in self.ag_structure_u:
                for chain in model:
                    if chain.get_id() in [self.heavy_chain_id_b,
                                          self.light_chain_id_b]:
                        model.detach_child(chain.get_id())

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

        for model in structure:
            for chain in model:
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
    def _load_structure(pdb_id, assembly_id):
        pdb = Conformation.pdb_parser.get_structure(pdb_id,
                                                    fetch_all_assemblies(
                                                        pdb_id)[
                                                        assembly_id - 1])

        tmp_path = os.path.join(DB_PATH, 'tmp.pdb')

        Conformation.pdb_io.set_structure(pdb)
        # delete all second variants from disordered atoms in order to get
        # rid of some problems
        Conformation.pdb_io.save(tmp_path, select=NotDisordered())

        pdb = Conformation.pdb_parser.get_structure(pdb_id, tmp_path)

        return union_models(pdb)

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
        # TODO: что-то не так тут выравнивается на 6icc
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
        path = os.path.join(DB_PATH, self.pdb_id_b, str(self.candidate_id))
        name_prefix = os.path.join(path,
                                   self.ab_pdb_id_u + '_' + self.ag_pdb_id_u)

        if not os.path.exists(path):
            os.makedirs(path)

        sb = StructureBuilder()

        sb.init_structure('complex')

        counter = 0

        for model in self.ab_structure_u.copy():
            model.id = counter
            sb.structure.add(model)
            counter += 1

        for model in self.ag_structure_u.copy():
            model.id = counter
            sb.structure.add(model)
            counter += 1

        models_in_struct = union_models(sb.structure)

        self.pdb_io.set_structure(models_in_struct)
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


def fetch_number_of_assemblies(pdb_id):
    curl = 'https://www.rcsb.org/pdb/rest/bioassembly/' \
           'nrbioassemblies?structureId={}' \
        .format(pdb_id)

    r = get_while_true(curl)
    xml = ElementTree.fromstring(r)

    return int(xml.attrib['count'])


def fetch_all_assemblies(pdb_id):
    n = fetch_number_of_assemblies(pdb_id)

    res = []

    for i in range(n):
        curl = 'https://files.rcsb.org/download/{}.pdb{}'. \
            format(pdb_id, str(i + 1))

        path_to_tmp = os.path.join(DB_PATH, pdb_id + '_BA_' + str(i) + DOT_PDB)

        if os.path.exists(path_to_tmp):
            res.append(path_to_tmp)
            continue

        r = get_while_true(curl)

        if r is None:
            continue

        with open(path_to_tmp, 'w') as f:
            f.write(r)

        res.append(path_to_tmp)

    return res


@memoize
def assembly_id_by_chains(pdb_id, chains):
    pdb_parser = PDBParser()

    counter = 1
    for assembly_path in fetch_all_assemblies(pdb_id):
        assembly_structure = pdb_parser.get_structure('ba', assembly_path)
        assembly = union_models(assembly_structure)

        chains_in_assembly = [x.get_id().split('_')[0]
                              for x in assembly.get_chains()]

        if frozenset(chains) <= frozenset(chains_in_assembly):
            return counter

        counter += 1

    return None


def union_models(struct):
    available_chain_ids = set(string.ascii_lowercase + string.ascii_uppercase)

    models = list(struct.get_models())

    if len(models) < 2:
        return struct

    sb = StructureBuilder()

    sb.init_structure('ba')
    sb.init_model(0)

    for model in models:
        for chain in model:
            chain_id = chain.get_id()
            name = chain_id

            if chain_id in available_chain_ids:
                available_chain_ids.remove(chain_id)
            else:
                # TODO: can crash if there are more than 52 chains
                name = available_chain_ids.pop()

            chain_copied = chain.copy()
            chain_copied.id = name

            sb.model.add(chain_copied)

    return sb.structure


class AssemblyMatchInfo:
    def __init__(self, assembly_id, matching, reason_bad=None):
        self.is_good = reason_bad is None
        self.id = assembly_id
        self.matching = matching
        self.reason_bad = reason_bad

    def __repr__(self):
        return str((self.is_good, self.id, self.matching,
                    self.reason_bad))


def check_structure(source_pdb_id, source_chain_ids, target_pdb_id, type):
    is_ab = type == AB

    source_seqs = list(filter(lambda x: x[0] in source_chain_ids,
                              fetch_all_sequences(source_pdb_id)))

    target_seqs = {k: v for k, v in fetch_all_sequences(target_pdb_id)}

    res = []

    n = -1
    for assembly_path in fetch_all_assemblies(target_pdb_id):
        n += 1

        pdb_parser = PDBParser()
        assembly_structure = pdb_parser.get_structure('ba', assembly_path)

        assembly = union_models(assembly_structure)

        chains_in_assembly = [x.get_id().split('_')[0]
                              for x in assembly.get_chains()]

        assembly_ids_seqs = list(
            map(lambda x: (x, target_seqs[x]), chains_in_assembly))

        chain_matching = defaultdict(list)

        for chain_id, chain_seq in source_seqs:
            for target_chain_id, target_seq in assembly_ids_seqs:
                if compare_query_and_hit_seqs(chain_seq, target_seq,
                                              None,
                                              None,
                                              write_log=False,
                                              is_ab=is_ab):
                    chain_matching[chain_id].append(target_chain_id)

        lens_of_matches = list(map(lambda x: len(chain_matching[x]),
                                   source_chain_ids))

        if not lens_of_matches:
            continue

        n_plus_one = n + 1

        if all(map(lambda x: x == 1, lens_of_matches)) and len(
                assembly_ids_seqs) == len(source_seqs):
            # assembly contains only matching with needed seqs

            res.append(AssemblyMatchInfo(n_plus_one, chain_matching))
        elif all(map(lambda x: x == 1, lens_of_matches)) and len(
                assembly_ids_seqs) != len(source_seqs):
            # assembly contains matching and some other chains

            res.append(AssemblyMatchInfo(n_plus_one, chain_matching,
                                         reason_bad='additional_chains'))
        elif lens_of_matches[0] > 0 and all(
                map(lambda x: x == lens_of_matches[0],
                    lens_of_matches)) and len(assembly_ids_seqs) == \
                lens_of_matches[0] * len(source_seqs):
            # assembly contains potential homomer that contains many matchings

            res.append(AssemblyMatchInfo(n_plus_one, chain_matching,
                                         reason_bad='potenial_homomer'))
        elif lens_of_matches[0] > 0 and all(
                map(lambda x: x == lens_of_matches[0],
                    lens_of_matches)) and len(assembly_ids_seqs) != \
                lens_of_matches[0] * len(source_seqs):
            # assembly contains potential complex homomer that
            # contains many matchings

            res.append(AssemblyMatchInfo(n_plus_one, chain_matching,
                                         reason_bad='potential_'
                                                    'complex_homomer'))

    return res


def get_pdb_ids(l, ty):
    return list(
        frozenset(map(lambda x: x[1], filter(lambda x: x[0] == ty, l))))


def matching_to_str(chains, matchings):
    n_matchings = len(matchings[list(matchings.keys())[0]])

    by_matching = []

    for i in range(n_matchings):
        by_matching.append([])
        for chain in chains:
            by_matching[-1].append(matchings[chain][i])

    return '|'.join(map(lambda x: ':'.join(x), by_matching))


def filter_candidates_pack(comp_name, pdb_id, candidate_pdb_ids, chain_ids, ty,
                           filtered_csv, rejected_csv):
    for candidate_pdb_id in candidate_pdb_ids:
        chains_str = ':'.join(chain_ids)

        assemblies = check_structure(pdb_id,
                                     chain_ids,
                                     candidate_pdb_id, ty)
        for assembly in assemblies:
            matching_str = matching_to_str(chain_ids, assembly.matching)

            if assembly.is_good:
                filtered_csv.write(','.join(
                    [comp_name, ty, chains_str, candidate_pdb_id,
                     matching_str, str(assembly.id)]) + '\n')
            else:
                rejected_csv.write(','.join(
                    [comp_name, ty, chains_str, candidate_pdb_id,
                     matching_str, str(assembly.id),
                     assembly.reason_bad]) + '\n')

        filtered_csv.flush()
        rejected_csv.flush()


def filter_for_unboundness(processed_csv):
    post_processed = set([])

    if os.path.exists('post_processed.csv'):
        with open('post_processed.csv', 'r') as f:
            for line in f.readlines():
                post_processed.add(line.strip())

    with open(FILTERED_COMPLEXES_CSV, 'w') as filtered_csv, open(
            REJECTED_COMPLEXES_CSV, 'w') as rejected_csv, open(
        'post_processed.csv', 'a') as post_processed_csv:
        filtered_csv.write(
            'comp_name,type,chain_ids,candidate_pdb_id,'
            'candidate_chain_ids,assembly_id\n')
        rejected_csv.write(
            'comp_name,type,chain_ids,candidate_pdb_id,'
            'candidate_chain_ids,assembly_id,reason\n')

        for comp_name, candidates in processed_csv.items():
            if comp_name in post_processed:
                continue

            pdb_id, chain_ids = comp_name.split('_')

            chain_ids = [x for x in chain_ids]
            heavy_chain_id = chain_ids[0]
            light_chain_id = chain_ids[1]
            ag_chains = chain_ids[2:]

            ab_candidates_pdb_ids = get_pdb_ids(candidates, AB)
            ag_candidates_pdb_ids = get_pdb_ids(candidates, AG)

            filter_candidates_pack(comp_name, pdb_id, ab_candidates_pdb_ids,
                                   [heavy_chain_id, light_chain_id], AB,
                                   filtered_csv, rejected_csv)

            filter_candidates_pack(comp_name, pdb_id, ag_candidates_pdb_ids,
                                   ag_chains, AG,
                                   filtered_csv, rejected_csv)

            post_processed_csv.write(str(comp_name) + '\n')
            post_processed_csv.flush()


def get_pbds_with_chains_and_assembly_ids(candidates, ty):
    return list(map(
        lambda x: (x.candidate_pdb_id, x.candidate_chain_ids, x.assembly_id),
        filter(lambda x: x.ty == ty, candidates)))


def process_candidates(comp_name, candidates):
    comp_name_split = comp_name.split('_')

    pdb_id_b = comp_name_split[0]
    chains_sep = list(map(lambda x: comp_name_split[1][x],
                          range(len(comp_name_split[1]))))

    # TODO: NO VHHs?
    heavy_chain_id_b = chains_sep[0]
    light_chain_id_b = chains_sep[1]
    ag_chain_ids_b = chains_sep[2:]

    ag_pdbs_with_chains = get_pbds_with_chains_and_assembly_ids(candidates, AG)
    ab_pdbs_with_chains = get_pbds_with_chains_and_assembly_ids(candidates, AB)

    is_ab_u = True
    is_ag_u = True

    assembly_id_b = assembly_id_by_chains(pdb_id_b,
                                          ag_chain_ids_b + [heavy_chain_id_b,
                                                            light_chain_id_b])

    # TODO: подумать тут
    if not assembly_id_b:
        return []

    if not ab_pdbs_with_chains:
        is_ab_u = False
        ab_pdbs_with_chains = [
            (pdb_id_b, heavy_chain_id_b + ':' + light_chain_id_b,
             assembly_id_b)]

    if not ag_pdbs_with_chains:
        is_ag_u = False
        ag_pdbs_with_chains = [
            (pdb_id_b, ':'.join(ag_chain_ids_b), assembly_id_b)]

    res = []

    counter = -1

    for ag_pdb_id_u, chains_ag, ag_assembly_id in ag_pdbs_with_chains:
        chains_ag_split = chains_ag.split(':')
        for ab_pdb_id_u, chains_ab, ab_assembly_id in ab_pdbs_with_chains:
            counter += 1

            [heavy_chain_id_u, light_chain_id_u] = chains_ab.split(':')
            conformation = Conformation(pdb_id_b, assembly_id_b,
                                        heavy_chain_id_b,
                                        light_chain_id_b,
                                        ag_chain_ids_b, ab_pdb_id_u,
                                        ab_assembly_id,
                                        heavy_chain_id_u,
                                        light_chain_id_u,
                                        ag_pdb_id_u, ag_assembly_id,
                                        chains_ag_split,
                                        is_ab_u, is_ag_u, counter)
            conformation.align_ab()
            conformation.align_ag()
            conformation.write_candidate()
            res.append(conformation)

    return res


def process_unbound(path_to_unbound_csv):
    prepared = process_csv(pd.read_csv(path_to_unbound_csv))

    for key, value in prepared.items():
        process_candidates(key, value)


class FilteredStructure:
    def __init__(self, line):
        self.ty = line['type']
        self.candidate_pdb_id = line['candidate_pdb_id']
        self.candidate_chain_ids = line['candidate_chain_ids']
        self.assembly_id = line['assembly_id']


def process_filtered_csv(path_to_filtered_structures_csv,
                         path_to_filtered_complexes_csv,
                         path_to_rejected_complexes_csv):
    filtered_structures_csv = pd.read_csv(path_to_filtered_structures_csv)

    by_complex = defaultdict(list)

    for i in range(len(filtered_structures_csv)):
        by_complex[filtered_structures_csv.iloc[i]['comp_name']].append(
            FilteredStructure(filtered_structures_csv.iloc[i]))

    with open(path_to_filtered_complexes_csv,
              'w') as filtered_complexes_csv, open(
        path_to_rejected_complexes_csv, 'w') as rejected_complexes_csv:
        for comp_name, structures in by_complex.items():
            process_candidates(comp_name, structures)


# process_unbound('unbound_data.csv')

# pdb_parser = PDBParser()
# assembly_structure = pdb_parser.get_structure('ba',
#                                               fetch_all_assemblies('1out')[0])
#
# assembly = union_models(assembly_structure)
#
# print(list(assembly.get_chains()))


# print(check_structure('6mfp', ['G'], '4dvv', AG))

# filter_for_unboundness(process_csv(pd.read_csv('unbound_data-3.csv')))
#
process_filtered_csv(FILTERED_STRUCTURES_CSV, FILTERED_COMPLEXES_CSV,
                     REJECTED_COMPLEXES_CSV)
