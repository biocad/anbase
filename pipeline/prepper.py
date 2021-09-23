import multiprocessing
import os
import random
import shutil
import subprocess
import time
import re
import pandas as pd

DB_PATH = 'data'
DOT_PDB = '.pdb'

HETATMS_DELETED = 'hetatms_deleted'

PDB_PREP_DIR = 'pdb_prep'
SCHROD_SCRIPT_PATH = 'resources/schrod_multi_prepare.sh'
PDBFIXER_SCRIPT_PATH = 'resources/pdbfixer_multi_prepare.sh'

PREP_SUFF = '.o.pdb'

PREPPED = 'prepared'

SCHROD = 'schrod'
PDBFIXER = 'PDBFixer'

SEQS = 'seqs'

DOT_FASTA = '.fasta'

DB_INFO_PATH = 'db_info_{}.csv'


def get_pdb_paths(dir_path, prev_epoch, comps_to_take, fast_version=False):
    def accept_path(p):
        if not comps_to_take:
            return True

        # m = re.search('/(...._(.\+.|.)\|(.\+.\+.\+.|.\+.\+.|.\+.|.))/', p)
        # comp_name = m.group(1)
        comp_name = os.path.basename(os.path.dirname(p))

        return comp_name in comps_to_take

    pdb_paths = []

    for root, dirnames, files in os.walk(dir_path):
        if not os.path.basename(root) == prev_epoch:
            continue

        if not accept_path(root):
            continue

        for file in files:
            if file.endswith(DOT_PDB):
                pdb_paths.append(os.path.join(root, file))

        random.shuffle(dirnames)
        candidate_dirs = dirnames if not fast_version else \
            dirnames[:min(2, len(dirnames))]

        for candidate_dir in candidate_dirs:
            candidate_dir_path = os.path.join(root, candidate_dir)

            for file in os.listdir(candidate_dir_path):
                if file.endswith(DOT_PDB):
                    pdb_paths.append(os.path.join(candidate_dir_path, file))

    return pdb_paths


def await_expected_files(expected_files, tmp_dir):
    res = {}

    waiting_time = 10 * 60

    while len(expected_files) > 0:
        for path in os.listdir('.'):
            if path.endswith(PREP_SUFF) and path in expected_files:
                expected_files.remove(path)
                name = path[:-len(PREP_SUFF)]
                res[name] = os.path.join(tmp_dir, path)
        time.sleep(1)
        waiting_time -= 1

        if waiting_time <= 0:
            break

    return res


def prep_in_mode(file_names, tmp_dir, mode, arg=None):
    expected_files = set(map(lambda x: x + PREP_SUFF, file_names))

    path_to_script = SCHROD_SCRIPT_PATH if mode == SCHROD else \
        PDBFIXER_SCRIPT_PATH
    script_name = os.path.basename(path_to_script)

    shutil.copyfile(path_to_script,
                    os.path.join(tmp_dir, script_name))
    shutil.copymode(path_to_script,
                    os.path.join(tmp_dir, script_name))

    os.chdir(tmp_dir)

    command = './{} {}'.format(script_name,
                               arg) if arg else './{}'.format(script_name)
    subprocess.call(command, stdout=subprocess.PIPE, shell=arg is not None)

    res = await_expected_files(expected_files, tmp_dir)

    os.chdir('..')

    return res


def schrod_prep(file_names, tmp_dir):
    return prep_in_mode(file_names, tmp_dir, SCHROD, DB_PATH)


def pdb_fixer_prep(file_names, tmp_dir):
    return prep_in_mode(file_names, tmp_dir, PDBFIXER, DB_PATH)


def prep_pdbs(last_epoch_name, epoch_name, db_path, mode, tmp_dir,
              fast_version, comps_to_take, run_id):
    pdbs_to_copy = []

    for file in os.listdir(db_path):
        dir_path = os.path.join(db_path, file)
        if os.path.isdir(dir_path):
            pdbs_to_copy += get_pdb_paths(dir_path, last_epoch_name, comps_to_take,
                                          fast_version=fast_version)

    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)

    name_to_path = {}

    for path in pdbs_to_copy:
        name = path.replace('+', '_').replace('|', '_').split('/')
        name.reverse()
        name = '_'.join(name).replace(DOT_PDB, '') + DOT_PDB

        name_to_path[name] = path

        shutil.copyfile(path.replace(last_epoch_name, SEQS).
                        replace(DOT_PDB, DOT_FASTA),
                        os.path.join(tmp_dir, name) + DOT_FASTA)
        shutil.copyfile(path, os.path.join(tmp_dir, name))

    unprepped_names = set(name_to_path.keys())

    while len(unprepped_names) > 0:
        print('New iteration:', len(unprepped_names), flush=True)

        with open('unprepped_{}_{}.log'.format(mode, run_id), 'w') as f:
            for name in unprepped_names:
                f.write(name_to_path[name] + '\n')
                f.flush()

        path_to_prepped_schrod = schrod_prep(unprepped_names,
                                             tmp_dir) if mode == SCHROD \
            else pdb_fixer_prep(unprepped_names, tmp_dir)

        for name in list(unprepped_names):
            if name not in path_to_prepped_schrod.keys():
                continue

            unprepped_names.remove(name)

            path_to_prepped = path_to_prepped_schrod[name]

            new_path = name_to_path[name].replace(last_epoch_name, epoch_name)

            if not os.path.exists(os.path.dirname(new_path)):
                os.makedirs(os.path.dirname(new_path))

            shutil.copyfile(path_to_prepped,
                            new_path)


if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('--db', default=DB_PATH, dest='db', metavar='DB',
                      help='Path to database [default: {}]'.format(DB_PATH))
    parser.add_option('--mode', default=SCHROD, dest='mode', metavar='MODE',
                      help='{} for preparation using schrodinger. '
                           '{} for preparation using PDBFixer. [default: {}]'.
                      format(SCHROD, PDBFIXER, SCHROD))
    parser.add_option('--prev-epoch', default=HETATMS_DELETED,
                      dest='prev_epoch', metavar='PREV_EPOCH',
                      help='Name of the epoch structures from which will be '
                           'prepared. [default: {}]'.format(HETATMS_DELETED))
    parser.add_option('--tmp-dir', default=PDB_PREP_DIR,
                      dest='tmp_dir', metavar='TMP_DIR',
                      help='Directory in which preparations will take place. '
                           'If it doesn\'t exist, it will be created '
                           'automatically.[default: {}]'.format(PDB_PREP_DIR))
    parser.add_option('--cur-epoch', default=PREPPED,
                      dest='cur_epoch', metavar='CUR_EPOCH',
                      help='Name of the preparation epoch. [default: {}]'.
                      format(PREPPED))
    parser.add_option('--fast-version', default=False,
                      dest='fast_version', metavar='FAST_VERSION',
                      help='Skip most of the alternative candidates. '
                           '[default: False]')
    parser.add_option('--run-id', default='0',
                      dest='run_id',
                      metavar='RUN_ID',
                      help='ID of the current run [default: {}]'.
                      format('0'))
    options, _ = parser.parse_args()

    if options.fast_version == 'True':
        is_fast_version = True
    else:
        is_fast_version = False

    comps_to_take = set()
    db_info = pd.read_csv(DB_INFO_PATH.format(options.run_id))

    for i in range(len(db_info)):
        comps_to_take.add(db_info.iloc[i]['comp_name'])

    prep_pdbs(options.prev_epoch, options.cur_epoch, options.db, options.mode,
              options.tmp_dir + '_' + options.run_id,
              is_fast_version, comps_to_take, options.run_id)
