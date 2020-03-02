import os
import pickle
from collections import defaultdict
from xml.etree import ElementTree

import numpy as np
import pandas as pd
from Bio import pairwise2
from Bio.PDB import PDBParser, Superimposer, PDBIO, Select
from Bio.PDB.Polypeptide import PPBuilder, is_aa
from Bio.PDB.StructureBuilder import StructureBuilder

from collect_db import AG, AB, DB_PATH, DOT_PDB, \
    fetch_sequence, memoize, get_while_true, \
    ANTIGEN_TYPE, PDB_ID, sub_nan, ANTIGEN_CHAIN, \
    H_CHAIN, L_CHAIN, form_comp_name, comp_name_to_pdb_and_chains, \
    fetch_all_sequences
from post_unboundness_filtering import union_models, rename_chains, \
    fetch_all_assemblies

FILTERED_STRUCTURES_CSV = 'filtered_for_unboundness.csv'
REJECTED_STRUCTURES_CSV = 'rejected_for_unboundness.csv'

FILTERED_COMPLEXES_CSV = 'filtered_complexes.csv'
REJECTED_COMPLEXES_CSV = 'rejected_complexes.csv'

ALIGNED_EPOCH = 'aligned'
HETATMS_DELETED = 'hetatms_deleted'

SEQUENCES = 'seqs'


def comp_name_to_dir_name(comp_name):
    return comp_name.replace(':', '+')


def dir_name_to_comp_name(dir_name):
    return dir_name.replace('+', ':')


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

    MAX_NUMBER_OF_ATOMS_IN_SM_TARGET = 5
    MAX_NUMBER_OF_ATOMS_IN_SM_COMMITMENT = 15

    def __init__(self, comp_name,
                 pdb_id_b, assembly_id_b, ab_chain_ids_b,
                 ag_chain_ids_b,
                 ab_pdb_id_u, ab_assembly_id, ab_chain_ids_u,
                 ag_pdb_id_u, ag_assembly_id, ag_chain_ids_u, is_ab_u, is_ag_u,
                 candidate_id):
        self.comp_name = comp_name

        self.pdb_id_b = pdb_id_b
        self.assembly_id_b = assembly_id_b
        self.ab_chain_ids_b = ab_chain_ids_b

        self.is_vhh = len(ab_chain_ids_b) == 1

        self.ag_chain_ids_b = ag_chain_ids_b

        self.ab_pdb_id_u = ab_pdb_id_u
        self.ab_chain_ids_u = ab_chain_ids_u

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
                chains = list(model.get_chains())
                for chain in chains:
                    if chain.get_id() not in self.ab_chain_ids_b:
                        model.detach_child(chain.get_id())

        if self.is_ag_u:
            self.ag_structure_u = self._load_structure(ag_pdb_id_u,
                                                       self.ag_assembly_id)
        else:
            self.ag_structure_u = self.complex_structure_b.copy()

            for model in self.ag_structure_u:
                chains = list(model.get_chains())
                for chain in chains:
                    if chain.get_id() in self.ab_chain_ids_b:
                        model.detach_child(chain.get_id())

        self.ab_chains_b = self.extract_chains(self.complex_structure_b,
                                               self.ab_chain_ids_b)
        self.ag_chains_b = self.extract_chains(self.complex_structure_b,
                                               self.ag_chain_ids_b)

        self.ab_atoms_b = []
        self.ag_atoms_b = []

        for chain in self.ab_chains_b:
            self.ab_atoms_b += self.extract_cas(chain)

        for chain in self.ag_chains_b:
            self.ag_atoms_b += self.extract_cas(chain)

        self.ab_interface_cas, self.ag_interface_cas = self.get_interface_cas()
        self.interface_atoms = list(self.ab_interface_cas) + list(
            self.ag_interface_cas)

        self.candidate_id = candidate_id

        self.is_aligned = False
        self.candidate_type = 'U:U' if self.is_ab_u and self.is_ag_u else \
            ('B:U' if self.is_ag_u else 'U:B')

        self.dir_name = comp_name_to_dir_name(self.comp_name)

        self.ab_seqs_b = None
        self.ag_seqs_b = None

        self.ab_seqs_u = None
        self.ag_seqs_u = None

    @staticmethod
    def extract_chains(structure, chain_ids):
        chains = []

        for chain_id in chain_ids:
            for model in structure:
                for chain in model:
                    if chain.get_id() == chain_id:
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

        return frozenset(ab_interface_cas), frozenset(ag_interface_cas)

    @staticmethod
    def _load_structure(pdb_id, assembly_id):
        assemblies = fetch_all_assemblies(pdb_id)
        pdb = Conformation.pdb_parser.get_structure(pdb_id,
                                                    assemblies[
                                                        assembly_id - 1])

        for x in assemblies:
            os.remove(x)

        tmp_path = os.path.join(DB_PATH, 'tmp.pdb')

        Conformation.pdb_io.set_structure(pdb)
        # delete all second variants from disordered atoms in order to get
        # rid of some problems
        Conformation.pdb_io.save(tmp_path, select=NotDisordered())

        pdb = Conformation.pdb_parser.get_structure(pdb_id, tmp_path)

        os.remove(tmp_path)

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

        mutual_ids = list(
            filter(lambda x: 'CA' in peps1[ids1[x]] and 'CA' in peps2[ids2[x]],
                   frozenset(ids1.keys()) & frozenset(ids2.keys())))

        atoms1 = [peps1[ids1[i]]['CA'] for i in mutual_ids]
        atoms2 = [peps2[ids2[i]]['CA'] for i in mutual_ids]

        return atoms1, atoms2

    @staticmethod
    def _inner_align(chain_ids_b, chains_b, pdb_id_b, structure_u, chain_ids_u,
                     pdb_id_u, atoms):
        chains_u = Conformation.extract_chains(structure_u, chain_ids_u)

        atoms1 = []
        atoms2 = []

        for i in range(len(chains_u)):
            tmp_atoms1, tmp_atoms2 = Conformation._matching_atoms_for_chains(
                chains_b[i],
                pdb_id_b,
                chain_ids_b[i],
                chains_u[i],
                pdb_id_u,
                chain_ids_u[i])

            atoms1 += tmp_atoms1
            atoms2 += tmp_atoms2

        interface_atoms1 = []
        interface_atoms2 = []

        for atom1, atom2 in zip(atoms1, atoms2):
            if atom1 in atoms:
                interface_atoms1.append(atom1)
                interface_atoms2.append(atom2)

        Conformation.super_imposer.set_atoms(interface_atoms1,
                                             interface_atoms2)
        Conformation.super_imposer.apply(structure_u.get_atoms())

        print(Conformation.super_imposer.rms)

    def _align_ab(self):
        if not self.is_ab_u:
            return

        self._inner_align(self.ab_chain_ids_b, self.ab_chains_b, self.pdb_id_b,
                          self.ab_structure_u, self.ab_chain_ids_u,
                          self.ab_pdb_id_u, self.ab_interface_cas)

    def _align_ag(self):
        if not self.is_ag_u:
            return

        self._inner_align(self.ag_chain_ids_b, self.ag_chains_b, self.pdb_id_b,
                          self.ag_structure_u, self.ag_chain_ids_u,
                          self.ag_pdb_id_u, self.ag_interface_cas)

    def alignment_epoch(self, epoch_name):
        self._align_ab()
        self._align_ag()
        self.write_candidate(epoch_name)

        self.is_aligned = True

    class SmallMoleculeStat:
        def __init__(self, name, n_atoms, dist):
            self.name = name
            self.n_atoms = n_atoms
            self.dist = dist

    def _get_small_molecule_stat_for_struct(self, small_molecules_stat_csv,
                                            molecule_res):

        atoms = list(molecule_res)

        mol_name = molecule_res.resname
        n_atoms = len(atoms)
        min_dist_to_interface = float('inf')

        for atom in atoms:
            for interface_atom in self.interface_atoms:
                min_dist_to_interface = min(
                    np.linalg.norm(atom.coord - interface_atom.coord),
                    min_dist_to_interface)

        if small_molecules_stat_csv:
            small_molecules_stat_csv.write('{},{},{},{},{},{:.2f}\n'.
                                           format(self.comp_name,
                                                  self.candidate_id,
                                                  self.candidate_type,
                                                  mol_name, n_atoms,
                                                  min_dist_to_interface))
            small_molecules_stat_csv.flush()

        return self.SmallMoleculeStat(mol_name, n_atoms, min_dist_to_interface)

    @staticmethod
    def _is_hoh(residue):
        return residue.resname == 'HOH'

    def _get_small_molecules_stat_for_struct(self, small_molecules_stat_csv,
                                             struct):
        non_aa_residues = []

        for model in struct:
            for chain in model:
                for residue in chain:
                    if not is_aa(residue) and not self._is_hoh(
                            residue):
                        non_aa_residues.append(residue)

        res = []

        for residue in non_aa_residues:
            res.append(self._get_small_molecule_stat_for_struct(
                small_molecules_stat_csv, residue))

        return res

    @staticmethod
    def are_good_mols(mols, a, b):
        return all(map(lambda x: x.n_atoms <= a or x.dist > b, mols))

    def get_small_molecules_stat(self, small_molecules_stat_csv):
        if not self.is_aligned:
            raise RuntimeError(
                'Small molecules statistics can be calculated only on aligned '
                'structures.')

        mols_ag = self._get_small_molecules_stat_for_struct(
            small_molecules_stat_csv,
            self.ag_structure_u)
        mols_ab = self._get_small_molecules_stat_for_struct(
            small_molecules_stat_csv,
            self.ab_structure_u)

        if small_molecules_stat_csv and len(mols_ab) == 0 and len(
                mols_ag) == 0:
            small_molecules_stat_csv.write('{},{},{},,,\n'.
                                           format(self.comp_name,
                                                  self.candidate_id,
                                                  self.candidate_type))
            small_molecules_stat_csv.flush()

        return mols_ag + mols_ab

    def small_molecules_logging(self, small_molecules_csv):
        mols = self.get_small_molecules_stat(None)

        if all(map(
                lambda x: x.n_atoms <= self.MAX_NUMBER_OF_ATOMS_IN_SM_TARGET,
                mols)):
            return

        message = None

        if all(map(lambda x: x.n_atoms <=
                             self.MAX_NUMBER_OF_ATOMS_IN_SM_COMMITMENT, mols)):
            message = 'small molecules with ' + \
                      str(self.MAX_NUMBER_OF_ATOMS_IN_SM_TARGET) + \
                      ' < n_atoms <= ' + \
                      str(self.MAX_NUMBER_OF_ATOMS_IN_SM_COMMITMENT) + \
                      ' detected'
        else:
            message = 'small molecules with' + \
                      ' n_atoms > ' + \
                      str(self.MAX_NUMBER_OF_ATOMS_IN_SM_COMMITMENT) + \
                      ' detected'

        small_molecules_csv.write('{},{},{},{}\n'.format(self.comp_name,
                                                         self.candidate_id,
                                                         self.candidate_type,
                                                         message))
        small_molecules_csv.flush()

    @staticmethod
    def delete_hetatms(struct):
        for model in struct:
            for chain in model:
                residues = list(chain)
                for residue in residues:
                    tags = residue.get_full_id()

                    if tags[3][0] != " ":
                        chain.detach_child(residue.get_id())

    def hetatms_deletion_epoch(self, epoch_name):
        self.delete_hetatms(self.ab_structure_u)
        self.delete_hetatms(self.ag_structure_u)
        self.delete_hetatms(self.complex_structure_b)
        self.write_candidate(epoch_name)

    @staticmethod
    def _load_sequences_for_pdb_and_chain_ids(prefix, pdb_id, name, chain_ids):
        all_seqs = fetch_all_sequences(pdb_id)

        res = []

        for chain_id in chain_ids:
            for seq_id, seq in all_seqs:
                if seq_id == chain_id:
                    path = os.path.join(prefix, '{}_{}.fasta'.format(name,
                                                                     seq_id))

                    res.append(seq)

                    if os.path.exists(path):
                        continue

                    with open(path, 'w') as f:
                        f.write('>' + seq_id + '\n')
                        f.write(seq + '\n')

                    break

        return res

    def load_sequences(self):
        dir_path = os.path.join(DB_PATH, self.dir_name, SEQUENCES)

        if not os.path.exists(dir_path):
            os.mkdir(dir_path)

        ab_seqs_b = self._load_sequences_for_pdb_and_chain_ids(dir_path,
                                                               self.pdb_id_b, self.pdb_id_b + '_r_b',
                                                               self.ab_chain_ids_b)

        if not self.ab_seqs_b:
            self.ab_seqs_b = ab_seqs_b

        ag_seqs_b = self._load_sequences_for_pdb_and_chain_ids(dir_path,
                                                               self.pdb_id_b, self.pdb_id_b + '_l_b',
                                                               self.ag_chain_ids_b)

        if not self.ag_seqs_b:
            self.ag_seqs_b = ag_seqs_b

        candidate_path = os.path.join(dir_path,
                                      str(self.candidate_id))

        if not os.path.exists(candidate_path):
            os.mkdir(candidate_path)

        name = self.ab_pdb_id_u + '_' + self.ag_pdb_id_u

        ab_seqs_u = self._load_sequences_for_pdb_and_chain_ids(candidate_path,
                                                               self.ab_pdb_id_u, name + '_r_u',
                                                               self.ab_chain_ids_u)

        if not self.ab_seqs_u:
            self.ab_seqs_u = ab_seqs_u

        ag_seqs_u = self._load_sequences_for_pdb_and_chain_ids(candidate_path,
                                                               self.ag_pdb_id_u, name + '_l_u',
                                                               self.ag_chain_ids_u)

        if not self.ag_seqs_u:
            self.ag_seqs_u = ag_seqs_u

    def load_candidate(self, epoch_name):
        prefix = os.path.join(DB_PATH, self.dir_name,
                              str(self.candidate_id), epoch_name,
                              self.ab_pdb_id_u + '_' + self.ag_pdb_id_u)

        self.ab_pdb_id_u = \
            self.pdb_parser.get_structure('receptor',
                                          prefix + '_r' + (
                                              '_u' if self.is_ab_u else '_b')
                                          + DOT_PDB)

        self.ab_pdb_id_u = \
            self.pdb_parser.get_structure('receptor',
                                          prefix + '_l' + (
                                              '_u' if self.is_ag_u else '_b')
                                          + DOT_PDB)

    def write_candidate(self, epoch_name):
        pre_path = os.path.join(DB_PATH, self.dir_name, epoch_name)

        if not os.path.exists(pre_path):
            os.makedirs(pre_path)

        complex_b_path = os.path.join(pre_path, self.pdb_id_b + DOT_PDB)

        self.pdb_io.set_structure(self.complex_structure_b)
        self.pdb_io.save(complex_b_path)

        path = os.path.join(pre_path, str(self.candidate_id))

        if not os.path.exists(path):
            os.makedirs(path)

        name_prefix = os.path.join(path,
                                   self.ab_pdb_id_u + '_' + self.ag_pdb_id_u)

        self.pdb_io.set_structure(rename_chains(self.ab_structure_u))
        self.pdb_io.save(
            name_prefix + '_r' + ('_u' if self.is_ab_u else '_b')
            + DOT_PDB)

        self.pdb_io.set_structure(rename_chains(self.ag_structure_u))
        self.pdb_io.save(
            name_prefix + '_l' + ('_u' if self.is_ag_u else '_b')
            + DOT_PDB)


def process_csv(csv):
    data = defaultdict(list)

    for i in range(len(csv)):
        data[csv.iloc[i]['comp_name']].append((csv.iloc[i]['type'],
                                               csv.iloc[i]['candidate_pdb_id'],
                                               csv.iloc[i][
                                                   'candidate_chain_names']))

    return data


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

        os.remove(assembly_path)

        counter += 1

    return None


def get_pbds_with_chains_and_assembly_ids(candidates, ty):
    return list(map(
        lambda x: (x.candidate_pdb_id, x.candidate_chain_ids, x.assembly_id),
        filter(lambda x: x.ty == ty, candidates)))


def get_candidates(comp_name, candidates, cache=False):
    pdb_id_b, ab_chain_ids_b, ag_chain_ids_b = comp_name_to_pdb_and_chains(
        comp_name)

    ag_pdbs_with_chains = get_pbds_with_chains_and_assembly_ids(candidates, AG)
    ab_pdbs_with_chains = get_pbds_with_chains_and_assembly_ids(candidates, AB)

    is_ab_u = True
    is_ag_u = True

    try:
        assembly_id_b = assembly_id_by_chains(pdb_id_b,
                                              ag_chain_ids_b + ab_chain_ids_b)
    except Exception as e:
        print('Couldn\'t process complex\'s assemblies:', comp_name, 'reason:',
              e, flush=True)
        return []

    # TODO: подумать тут
    if not assembly_id_b:
        return []

    if not ab_pdbs_with_chains:
        is_ab_u = False
        ab_pdbs_with_chains = [
            (pdb_id_b, ':'.join(ab_chain_ids_b),
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

            ab_chain_ids_u = chains_ab.split(':')
            try:
                pickled_path = os.path.join(DB_PATH,
                                            comp_name.replace(':', '+') + '_'
                                            + str(counter) + '.pickle')

                conformation = None

                if cache and os.path.exists(pickled_path):
                    with open(pickled_path, 'rb') as f:
                        conformation = pickle.load(f)
                else:
                    conformation = Conformation(comp_name, pdb_id_b,
                                                assembly_id_b,
                                                ab_chain_ids_b,
                                                ag_chain_ids_b, ab_pdb_id_u,
                                                ab_assembly_id,
                                                ab_chain_ids_u,
                                                ag_pdb_id_u, ag_assembly_id,
                                                chains_ag_split,
                                                is_ab_u, is_ag_u, counter)

                    if cache:
                        with open(pickled_path, 'wb') as f:
                            pickle.dump(conformation, f)

                res.append(conformation)
            except Exception as e:
                print('Couldn\'t process candidate', pdb_id_b, ab_pdb_id_u,
                      'assembly', ab_assembly_id, ag_pdb_id_u, 'assembly',
                      ag_assembly_id, e, flush=True)

    return res


class FilteredStructure:
    def __init__(self, line):
        self.ty = line['type']
        self.candidate_pdb_id = line['candidate_pdb_id']
        self.candidate_chain_ids = line['candidate_chain_ids']
        self.assembly_id = line['assembly_id']


def process_filtered_csv(path_to_filtered_structures_csv,
                         path_to_rejected_complexes_csv, to_accept=None):
    filtered_structures_csv = pd.read_csv(path_to_filtered_structures_csv)

    by_complex = defaultdict(list)

    for i in range(len(filtered_structures_csv)):
        by_complex[filtered_structures_csv.iloc[i]['comp_name']].append(
            FilteredStructure(filtered_structures_csv.iloc[i]))

    filter_out_peptides(by_complex,
                        pd.read_csv('data/sabdab_summary_all.tsv', sep='\t'))

    with open(path_to_rejected_complexes_csv, 'w') as rejected_complexes_csv, \
            open('small_molecules_log.csv', 'w') as small_molecules_log_csv:

        small_molecules_log_csv.write(
            'comp_name,candidate_id,candidate_type,reason\n')
        small_molecules_log_csv.flush()

        rejected_complexes_csv.write(
            'comp_name,candidate_id,candidate_type,reason\n')
        rejected_complexes_csv.flush()

        with_candidates = {}

        for comp_name, structures in list(by_complex.items()):
            if to_accept and comp_name not in to_accept:
                continue

            with_candidates[comp_name] = get_candidates(comp_name, structures,
                                                        cache=True)

        counter = 1
        for comp_name, candidates in with_candidates.items():
            print('Processing complex', comp_name,
                  '[{}/{}]'.format(counter, len(with_candidates)))

            for candidate in candidates:
                try:
                    candidate.alignment_epoch(ALIGNED_EPOCH)

                    candidate.small_molecules_logging(small_molecules_log_csv)

                    candidate.hetatms_deletion_epoch(HETATMS_DELETED)

                    candidate.load_sequences()

                except Exception as e:
                    rejected_complexes_csv.write('{},{},{},{}\n'.format(
                        candidate.comp_name, candidate.candidate_id,
                        candidate.candidate_type, str(e)))
                    rejected_complexes_csv.flush()

            counter += 1


def filter_out_peptides(filtered_structures, sabdab_tb):
    peptide_complexes = set([])

    for i in range(len(sabdab_tb)):
        entry = sabdab_tb.iloc[i]
        antigen_type = sub_nan(entry[ANTIGEN_TYPE])

        if antigen_type and 'peptide' in antigen_type:
            comp_name = form_comp_name(entry[PDB_ID],
                                       [sub_nan(entry[H_CHAIN]),
                                        sub_nan(entry[L_CHAIN])], entry[
                                           ANTIGEN_CHAIN].split(' | '))
            peptide_complexes.add(comp_name)

    for x in peptide_complexes:
        if x in filtered_structures:
            del filtered_structures[x]


if __name__ == '__main__':
    process_filtered_csv(FILTERED_STRUCTURES_CSV,
                         REJECTED_COMPLEXES_CSV)
