import os
from collections import defaultdict
from multiprocessing.pool import Pool

from Bio import pairwise2
from Bio.PDB import PDBParser

from candidate_info import CandidateInfo
from fetch_unbound_data import comp_name_to_pdb_and_chains
from process_unbound_data import SEQUENCES, HETATMS_DELETED, \
    INTERFACE_CUTOFF, extract_seq
import pandas as pd
import numpy as np

from candidate_info import read_fasta
from fetch_unbound_data import fetch_all_sequences
from process_unbound_data import Conformation

DB_PATH = 'data'
DOT_PDB = '.pdb'
DOT_FASTA = '.fasta'

DB_INFO_PATH = 'db_info.csv'

pdb_parser = PDBParser()

CONSTRAINTS = 'constraints'
EPITOPE = 'epitope'

CLOSE_CUTOFF = 6.5


def generate_constraints_task(conformation_like, epoch_name):
    candidate_name = conformation_like.comp_name + '_' + \
                     str(conformation_like.candidate_id)
    print('Processing candidate:', candidate_name)

    try:
        generate_constraints(conformation_like, epoch_name)
    except Exception as e:
        print('Couldn\'t process candidate:', candidate_name,
              'reason:', e, flush=True)


def generate_constraints(conformation_like, epoch_name):
    def group(l):
        res = []

        for x in l:
            if len(res) == 0:
                res.append((x, x))
            elif x.isnumeric() and res[-1][1].isnumeric() \
                    and int(x) == int(res[-1][1]) + 1:
                res[-1] = (res[-1][0], x)
            else:
                res.append((x, x))

        return res

    def form_constraints(atoms_set, chain_ids):
        chains_to_constraints = {x: [] for x in chain_ids}

        atoms = list(atoms_set)

        atoms.sort(key=lambda x: x.get_full_id())

        for atom in atoms:
            _, residue_id, residue_suf = atom.get_parent().get_id()

            chain_id = atom.get_parent().get_parent().get_id()
            res_name = str(residue_id) + residue_suf.strip()

            if len(chains_to_constraints[chain_id]) == 0:
                chains_to_constraints[chain_id].append(res_name)
            elif chains_to_constraints[chain_id][-1] != res_name:
                chains_to_constraints[chain_id].append(res_name)

        for key in chains_to_constraints:
            chains_to_constraints[key] = group(
                chains_to_constraints[key])

        return chains_to_constraints

    def write_constraints(path, constraints):
        with open(path, 'w') as f:
            for chain_id, ranges in constraints.items():
                f.write('>{}:{}\n'.format(chain_id, 'attraction'))
                f.write(','.join(map(lambda x: x[0] if x[0] == x[1] else
                '{}-{}'.format(x[0], x[1]), ranges)) + '\n')

    _, ag_interface_atoms_b = Conformation.get_interface_atoms(
        conformation_like.comp_name,
        conformation_like.ab_chains_b,
        conformation_like.ag_chains_b,
        dist=CLOSE_CUTOFF,
        only_ca=False)
    _, ag_interface_atoms_u = Conformation.get_corresponding_atoms(
        conformation_like.ag_chain_ids_b, conformation_like.ag_chains_b,
        conformation_like.pdb_id_b,
        conformation_like.ag_structure_u, conformation_like.ag_chain_ids_u,
        conformation_like.ag_pdb_id_u, ag_interface_atoms_b,
        only_cas=False,
        seqs_b=conformation_like.ag_seqs_b,
        seqs_u=conformation_like.ag_seqs_u)

    pre_path = os.path.join(DB_PATH, conformation_like.comp_name, epoch_name)

    chains_to_constraints_b = form_constraints(ag_interface_atoms_b,
                                               conformation_like.ag_chain_ids_b)
    chains_to_constraints_u = form_constraints(ag_interface_atoms_u,
                                               conformation_like.ag_chain_ids_u)

    path_to_constraints_b = os.path.join(pre_path,
                                         conformation_like.pdb_id_b +
                                         '_ag_b.fasta')

    if not os.path.exists(os.path.dirname(path_to_constraints_b)):
        os.makedirs(os.path.dirname(path_to_constraints_b))

    write_constraints(path_to_constraints_b, chains_to_constraints_b)

    path_to_candidate = os.path.join(pre_path,
                                     str(conformation_like.candidate_id),
                                     conformation_like.pdb_id_b +
                                     '_ag_u.fasta')

    if not os.path.exists(os.path.dirname(path_to_candidate)):
        os.makedirs(os.path.dirname(path_to_candidate))

    write_constraints(path_to_candidate, chains_to_constraints_u)


if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('--db', default=DB_PATH, dest='db', metavar='DB',
                      help='Path to database [default: {}]'.format(DB_PATH))
    parser.add_option('--db-info', default=DB_INFO_PATH, dest='db_info',
                      metavar='DB_INFO_PATH',
                      help='Path to database info csv file [default: {}]'.
                      format(DB_INFO_PATH))
    parser.add_option('--prev-epoch', default=HETATMS_DELETED,
                      dest='prev_epoch', metavar='PREV_EPOCH',
                      help='Name of the epoch structures from which will be '
                           'checked for gaps. [default: {}]'.format(
                          HETATMS_DELETED))
    parser.add_option('--name-of-constraints-folder', default=EPITOPE,
                      dest='constraints_folder_name',
                      metavar='CONSTRAINTS_FOLDER_NAME',
                      help='Name of the folder where constraints will be '
                           'stored for each complex. [default: {}]'.format(
                          EPITOPE))
    parser.add_option('--only-uu', default=False,
                      dest='only_uu', metavar='ONLY_UU',
                      help='Flag to process only candidates of type UU. '
                           '[default: False]')
    parser.add_option('--number-of-processes', default=1,
                      dest='number_of_processes',
                      metavar='NUMBER_OF_PROCESSES',
                      help='Number of processes to use. '
                           '[default: 1]')
    options, _ = parser.parse_args()

    df = pd.read_csv(options.db_info, dtype=str)

    with Pool(int(options.number_of_processes)) as pool:
        tasks = []

        for i in range(len(df)):
            candidate_info = CandidateInfo(df.iloc[i])

            print('Reading', candidate_info.comp_name + '_' + str(
                candidate_info.candidate_id), '[{}/{}]'.format(i + 1, len(df)),
                  flush=True)

            if options.only_uu and candidate_info.candidate_type != 'U:U':
                continue

            try:
                if int(options.number_of_processes) == 1:
                    generate_constraints(
                        candidate_info.to_conformation_like(options.db,
                                                            options.prev_epoch),
                        options.constraints_folder_name)
                else:
                    tasks.append(
                        (candidate_info.to_conformation_like(options.db,
                                                             options.prev_epoch),
                         options.constraints_folder_name))
            except Exception as e:
                print(e)

        if len(tasks) > 0:
            pool.starmap(generate_constraints_task, tasks)
