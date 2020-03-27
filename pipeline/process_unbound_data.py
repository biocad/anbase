import os
import pickle
import string
from collections import defaultdict
from xml.etree import ElementTree

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, Superimposer, PDBIO, Select
from Bio.PDB.Polypeptide import PPBuilder, is_aa, d1_to_index, dindex_to_3

import alignments
from fetch_unbound_data import AG, AB, DB_PATH, DOT_PDB, \
    fetch_sequence, memoize, \
    ANTIGEN_TYPE, PDB_ID, sub_nan, ANTIGEN_CHAIN, \
    H_CHAIN, L_CHAIN, form_comp_name, comp_name_to_pdb_and_chains, \
    fetch_all_sequences, CHAINS_SEPARATOR, \
    calc_mismatches_stat, extract_seq, retrieve_resolution
from filter_unbound_data import union_models, \
    fetch_all_assemblies

FILTERED_STRUCTURES_CSV = 'filtered_for_unboundness.csv'
REJECTED_STRUCTURES_CSV = 'rejected_for_unboundness.csv'

FILTERED_COMPLEXES_CSV = 'filtered_complexes.csv'
REJECTED_COMPLEXES_CSV = 'rejected_complexes.csv'

ALIGNED = 'aligned'
HETATMS_DELETED = 'hetatms_deleted'

SEQUENCES = 'seqs'

INTERFACE_CUTOFF = 10

DB_INFO_COLUMNS = ['comp_name', 'candidate_type', 'candidate_id',
                   'pdb_id_b', 'ab_chain_ids_b', 'ag_chain_ids_b',
                   'resolution_b', 'resolution_method_b',
                   'ab_pdb_id_u', 'ab_chain_ids_u', 'ab_resolution_u',
                   'ab_resolution_method_u',
                   'ag_pdb_id_u', 'ag_chain_ids_u', 'ag_resolution_u',
                   'ag_resolution_method_u',
                   'ab_mismatches_cnt', 'ag_mismatches_cnt',
                   'small_molecules_message']
DB_INFO_HEADER = ','.join(DB_INFO_COLUMNS)

REJECTED_COMPLEXES_COLUMNS = ['comp_name', 'candidate_id', 'candidate_type',
                              'reason']
REJECTED_COMPLEXES_HEADER = ','.join(REJECTED_COMPLEXES_COLUMNS)


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

    MAX_NUMBER_OF_ATOMS_IN_SM_TARGET = 7
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

        self.ab_assembly_id_u = ab_assembly_id
        self.ag_assembly_id_u = ag_assembly_id

        self.is_ab_u = is_ab_u

        self.is_ag_u = is_ag_u

        self.candidate_id = candidate_id

        self.complex_structure_b = self._load_structure(pdb_id_b,
                                                        self.assembly_id_b)
        if self.is_ab_u:
            self.ab_structure_u = self._load_structure(ab_pdb_id_u,
                                                       self.ab_assembly_id_u)
        else:
            self.ab_structure_u = self.complex_structure_b.copy()

            for model in self.ab_structure_u:
                chains = list(model.get_chains())
                for chain in chains:
                    if chain.get_id() not in self.ab_chain_ids_b:
                        model.detach_child(chain.get_id())

        if self.is_ag_u:
            self.ag_structure_u = self._load_structure(ag_pdb_id_u,
                                                       self.ag_assembly_id_u)
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

        self.complex_structure_b, self.complex_mapping_b = rename_chains(
            self.complex_structure_b)
        self.ab_structure_u, self.ab_mapping_u = rename_chains(
            self.ab_structure_u)
        self.ag_structure_u, self.ag_mapping_u = rename_chains(
            self.ag_structure_u)

        self.ab_atoms_b = []
        self.ag_atoms_b = []

        for chain in self.ab_chains_b:
            self.ab_atoms_b += self.extract_cas(chain)

        for chain in self.ag_chains_b:
            self.ag_atoms_b += self.extract_cas(chain)

        self.ab_interface_cas, self.ag_interface_cas = self.get_interface_cas()
        self.interface_atoms = list(self.ab_interface_cas) + list(
            self.ag_interface_cas)

        self.is_aligned = False
        self.candidate_type = 'U:U' if self.is_ab_u and self.is_ag_u else \
            ('B:U' if self.is_ag_u else 'U:B')

        self.dir_name = self.comp_name

        self.ab_seqs_b = None
        self.ag_seqs_b = None

        self.ab_seqs_u = None
        self.ag_seqs_u = None

    @staticmethod
    def prepend_sequence_info_to_pdb(pdb_path, pdb_id, mapping):
        all_seqs = fetch_all_sequences(pdb_id)

        def up_to(i, n):
            res = str(i)

            while len(res) < n:
                res = ' ' + res

            return res

        def to_3(x, i, seq_len):
            if x == 'X' and i == 0:
                return 'ACE'
            elif x == 'X' and i == seq_len - 1:
                return 'NME'
            elif x in d1_to_index:
                return dindex_to_3[d1_to_index[x]]
            else:
                print('WHAT AA is this:', pdb_path, pdb_id, i, flush=True)
                return 'UNK'

        def seq_to_seqres_section(seq, chain_name):
            seqres = 'SEQRES'
            n_of_residue_columns = 13
            len_of_pdb_row = 80

            full_names = list(
                map(lambda x: to_3(x[0], x[1], len(seq)),
                    zip(list(seq), range(len(seq)))))

            rows = []

            i = 0
            while len(full_names) > 0:
                i += 1

                to_take = min(n_of_residue_columns, len(full_names))

                row_names = full_names[:to_take]
                full_names = full_names[to_take:]

                ser_num = up_to(i, 3)
                num_res = up_to(len(seq), 4)

                row = '{} {} {} {}  {}'.format(seqres, ser_num,
                                               chain_name,
                                               num_res, ' '.join(row_names))

                rows.append(row + (len_of_pdb_row - len(row)) * ' ')

            return '\n'.join(rows)

        seqres_info = '\n'.join(
            map(lambda x: seq_to_seqres_section(all_seqs[mapping[x]], x),
                mapping.keys()))

        with open(pdb_path, 'r') as f:
            lines = f.readlines()

        with open(pdb_path, 'w') as f:
            f.write(seqres_info + '\n')
            f.writelines(lines)

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
        ab_interface_cas = []
        ag_interface_cas = []

        for ab_at in self.ab_atoms_b:
            for ag_at in self.ag_atoms_b:
                if np.linalg.norm(
                        ab_at.coord - ag_at.coord) < INTERFACE_CUTOFF:
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
        seq1 = fetch_sequence(pdb_id1, chain_id1)
        seq2 = fetch_sequence(pdb_id2, chain_id2)

        return Conformation._matching_atoms_for_chains_seqs(chain1, seq1,
                                                            chain2, seq2)

    @staticmethod
    def _matching_atoms_for_chains_seqs(chain1, seq1, chain2, seq2):
        def extract_peps(chain):
            peps = []

            for x in Conformation.peptides_builder.build_peptides(chain):
                peps += x

            return peps

        def get_local_ids_from_chain(chain, seq, ids_in_seq):
            struct_seq = extract_seq(chain)

            alignment_loc = \
                alignments.align_possibly_gapped_sequence_on_its_complete_version(
                    struct_seq, seq)[0]

            counter_local = -1
            counter_seq = -1

            res = []

            for i in range(len(alignment_loc[0])):
                if alignment_loc[0][i] == '-' and alignment_loc[1][i] == '-':
                    continue
                elif alignment_loc[0][i] == '-':
                    counter_seq += 1
                    continue
                elif alignment_loc[1][i] == '-':
                    counter_local += 1
                    continue
                else:
                    counter_local += 1
                    counter_seq += 1

                if counter_seq in ids_in_seq:
                    res.append((counter_seq, counter_local))

            return {key: value for key, value in res}

        alignment = \
            alignments.subsequence_without_gaps(
                seq1, seq2)[0]

        counter1 = -1
        counter2 = -1

        seq1_to_universal = {}
        seq2_to_universal = {}

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

            seq1_to_universal[counter1] = i
            seq2_to_universal[counter2] = i

        peps1 = extract_peps(chain1)
        peps2 = extract_peps(chain2)

        universal_to_seq1 = {k: v for v, k in seq1_to_universal.items()}
        universal_to_seq2 = {k: v for v, k in seq2_to_universal.items()}

        seq1_to_local = get_local_ids_from_chain(chain1, seq1,
                                                 seq1_to_universal.keys())
        seq2_to_local = get_local_ids_from_chain(chain2, seq2,
                                                 seq2_to_universal.keys())

        common_universal_ids = frozenset(map(lambda x: seq1_to_universal[x],
                                             seq1_to_local.keys())) & \
                               frozenset(map(lambda x: seq2_to_universal[x],
                                             seq2_to_local.keys()))

        ca_common_universal_ids = list(
            filter(
                lambda x: 'CA' in peps1[
                    seq1_to_local[universal_to_seq1[x]]] and 'CA' in
                          peps2[seq2_to_local[universal_to_seq2[x]]],
                common_universal_ids))

        atoms1 = [peps1[seq1_to_local[universal_to_seq1[i]]]['CA'] for i in
                  ca_common_universal_ids]
        atoms2 = [peps2[seq2_to_local[universal_to_seq2[i]]]['CA'] for i in
                  ca_common_universal_ids]

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

    def _load_sequences_for_pdb_and_chain_ids(self, prefix, pdb_id, name,
                                              chain_ids,
                                              mapping):
        path = os.path.join(prefix, '{}.fasta'.format(name))

        if os.path.exists(path):
            res = [None for _ in chain_ids]

            mapping = {k: v for k, v in zip(chain_ids, range(len(chain_ids)))}

            with open(path, 'r') as f:
                lines = list(map(lambda x: x.strip(), f.readlines()))

                i = 0

                while i < len(lines):
                    chain_id = lines[i].split(':')[1]
                    chain_seq = lines[i + 1]

                    if chain_id in chain_ids:
                        res[mapping[chain_id]] = chain_seq
                    i += 2

            return res

        all_seqs = fetch_all_sequences(pdb_id)

        res = []

        for new_chain_name, old_chain_name in mapping.items():
            with open(path, 'a') as f:
                f.write(
                    '>{}:{}\n'.format(self.pdb_id_b.upper(), new_chain_name))
                f.write(all_seqs[old_chain_name] + '\n')

        for chain_id in chain_ids:
            for seq_id, seq in all_seqs.items():
                if seq_id == chain_id:
                    res.append(seq)
                    break

        return res

    def load_sequences(self):
        dir_path = os.path.join(DB_PATH, self.dir_name, SEQUENCES)

        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

        seqs_b = self. \
            _load_sequences_for_pdb_and_chain_ids(dir_path,
                                                  self.pdb_id_b,
                                                  self.pdb_id_b,
                                                  self.ab_chain_ids_b +
                                                  self.ag_chain_ids_b,
                                                  self.complex_mapping_b)

        if not self.ab_seqs_b:
            self.ab_seqs_b = seqs_b[:len(self.ab_chain_ids_b)]

        if not self.ag_seqs_b:
            self.ag_seqs_b = seqs_b[len(self.ab_chain_ids_b):]

        candidate_path = os.path.join(dir_path,
                                      str(self.candidate_id))

        if not os.path.exists(candidate_path):
            os.mkdir(candidate_path)

        name = self.pdb_id_b

        ab_seqs_u = self._load_sequences_for_pdb_and_chain_ids(candidate_path,
                                                               self.ab_pdb_id_u,
                                                               name + '_ab_u',
                                                               self.ab_chain_ids_u,
                                                               self.ab_mapping_u)

        if not self.ab_seqs_u:
            self.ab_seqs_u = ab_seqs_u

        ag_seqs_u = self._load_sequences_for_pdb_and_chain_ids(candidate_path,
                                                               self.ag_pdb_id_u,
                                                               name + '_ag_u',
                                                               self.ag_chain_ids_u,
                                                               self.ag_mapping_u)

        if not self.ag_seqs_u:
            self.ag_seqs_u = ag_seqs_u

    @staticmethod
    def write_structure(structure, path, pdb_id, mapping):
        Conformation.pdb_io.set_structure(structure)
        Conformation.pdb_io.save(path)
        Conformation.prepend_sequence_info_to_pdb(path, pdb_id, mapping)

    @staticmethod
    def filter_chains(struct, chain_ids):
        for model in struct:
            chains = list(model.get_chains())
            for chain in chains:
                if chain.get_id() not in chain_ids:
                    model.detach_child(chain.get_id())

    def write_candidate(self, epoch_name):
        pre_path = os.path.join(DB_PATH, self.dir_name, epoch_name)

        if not os.path.exists(pre_path):
            os.makedirs(pre_path)

        complex_b_path = os.path.join(pre_path, self.pdb_id_b + DOT_PDB)

        complex_ab_b_path = os.path.join(pre_path,
                                         self.pdb_id_b + '_ab_b' + DOT_PDB)
        complex_ag_b_path = os.path.join(pre_path,
                                         self.pdb_id_b + '_ag_b' + DOT_PDB)

        if not os.path.exists(complex_b_path):
            self.write_structure(self.complex_structure_b, complex_b_path,
                                 self.pdb_id_b, self.complex_mapping_b)

        if not os.path.exists(complex_ab_b_path):
            ab_struct_b = self.complex_structure_b.copy()
            self.filter_chains(ab_struct_b, self.ab_chain_ids_b)

            self.write_structure(ab_struct_b, complex_ab_b_path,
                                 self.pdb_id_b,
                                 {x: x for x in self.ab_chain_ids_b})

        if not os.path.exists(complex_ag_b_path):
            ag_struct_b = self.complex_structure_b.copy()
            self.filter_chains(ag_struct_b, self.ag_chain_ids_b)

            self.write_structure(ag_struct_b, complex_ag_b_path,
                                 self.pdb_id_b,
                                 {x: x for x in self.ag_chain_ids_b})

        path = os.path.join(pre_path, str(self.candidate_id))

        if not os.path.exists(path):
            os.makedirs(path)

        name_prefix = os.path.join(path, self.pdb_id_b)

        self.write_structure(self.ab_structure_u, name_prefix + '_ab_u'
                             + DOT_PDB, self.ab_pdb_id_u, self.ab_mapping_u)

        self.write_structure(self.ag_structure_u, name_prefix + '_ag_u'
                             + DOT_PDB, self.ag_pdb_id_u, self.ag_mapping_u)

    MOLS_WARNING = 'small molecules with ' + \
                   str(MAX_NUMBER_OF_ATOMS_IN_SM_TARGET) + \
                   ' < n_atoms <= ' + \
                   str(MAX_NUMBER_OF_ATOMS_IN_SM_COMMITMENT) + \
                   ' detected'

    MOLS_ERROR = 'small molecules with' + \
                 ' n_atoms > ' + \
                 str(MAX_NUMBER_OF_ATOMS_IN_SM_COMMITMENT) + \
                 ' detected'

    @staticmethod
    def calc_mismatches(seqs1, seqs2):
        cnt = 0

        for x, y in zip(seqs1, seqs2):
            _, c, _ = calc_mismatches_stat(x, y)
            cnt += c

        return cnt

    def write_info(self, db_info_csv):
        mols = self.get_small_molecules_stat(None)

        if all(map(
                lambda x: x.n_atoms <= self.MAX_NUMBER_OF_ATOMS_IN_SM_TARGET,
                mols)):
            mols_message = 'NA'
        elif all(map(lambda x: x.n_atoms <=
                               self.MAX_NUMBER_OF_ATOMS_IN_SM_COMMITMENT,
                     mols)):
            mols_message = self.MOLS_WARNING
        else:
            mols_message = self.MOLS_ERROR

        resolution_b, resolution_method_b = retrieve_resolution(self.pdb_id_b)
        ab_resolution_u, ab_resolution_method_u = retrieve_resolution(
            self.ab_pdb_id_u)
        ag_resolution_u, ag_resolution_method_u = retrieve_resolution(
            self.ag_pdb_id_u)

        ab_mismatches = self.calc_mismatches(self.ab_seqs_u, self.ab_seqs_b)
        ag_mismatches = self.calc_mismatches(self.ag_seqs_u, self.ag_seqs_b)

        db_info_csv.write(','.join(['{}'] * len(DB_INFO_COLUMNS)).format(
            self.comp_name, self.candidate_type, self.candidate_id,
            self.pdb_id_b, CHAINS_SEPARATOR.join(self.ab_chain_ids_b),
            CHAINS_SEPARATOR.join(self.ag_chain_ids_b), resolution_b,
            resolution_method_b,
            self.ab_pdb_id_u, CHAINS_SEPARATOR.join(self.ab_chain_ids_u),
            ab_resolution_u, ab_resolution_method_u,
            self.ag_pdb_id_u, CHAINS_SEPARATOR.join(self.ag_chain_ids_u),
            ag_resolution_u, ag_resolution_method_u,
            ab_mismatches, ag_mismatches, mols_message) + '\n')
        db_info_csv.flush()


def rename_chains(struct):
    new_struct = struct.copy()

    available_chain_ids = set(
        string.ascii_lowercase + string.ascii_uppercase + string.digits)

    mapping = {}

    for model in new_struct:
        for chain in model:
            chain_id = chain.get_id()

            if chain_id in available_chain_ids:
                mapping[chain_id] = chain_id
                available_chain_ids.remove(chain_id)
            else:
                chain.id = available_chain_ids.pop()

                mapping[chain.id] = chain_id[0]

    return new_struct, mapping


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
            (pdb_id_b, CHAINS_SEPARATOR.join(ab_chain_ids_b),
             assembly_id_b)]

    if not ag_pdbs_with_chains:
        is_ag_u = False
        ag_pdbs_with_chains = [
            (pdb_id_b, CHAINS_SEPARATOR.join(ag_chain_ids_b), assembly_id_b)]

    res = []

    counter = -1

    for ag_pdb_id_u, chains_ag, ag_assembly_id in ag_pdbs_with_chains:
        chains_ag_split = chains_ag.split(CHAINS_SEPARATOR)
        for ab_pdb_id_u, chains_ab, ab_assembly_id in ab_pdbs_with_chains:
            counter += 1

            ab_chain_ids_u = chains_ab.split(CHAINS_SEPARATOR)
            try:
                pickled_path = os.path.join(DB_PATH,
                                            comp_name + '_'
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
            open('db_info.csv', 'w') as db_info_csv:

        db_info_csv.write(DB_INFO_HEADER + '\n')
        db_info_csv.flush()

        rejected_complexes_csv.write(REJECTED_COMPLEXES_HEADER + '\n')
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
                    candidate.load_sequences()
                    candidate.alignment_epoch(ALIGNED)
                    candidate.write_info(db_info_csv)
                    candidate.hetatms_deletion_epoch(HETATMS_DELETED)

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
