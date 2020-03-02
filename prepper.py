import os
import shutil
import subprocess
import time

from collect_db import DB_PATH, DOT_PDB
from post_hetatms_and_mols_elimination import HETATMS_DELETED

PDB_PREP_DIR = 'pdb_prep'
SCHROD_SCRIPT_PATH = 'schrod_multi_prepare.sh'

SCHROD_PREP_SUFF = '.o.pdb'

PREPPED = 'prepared'


def get_pdb_paths(dir_path, prev_epoch):
    pdb_paths = []

    for root, _, files in os.walk(dir_path):
        for file in files:
            if prev_epoch in root and file.endswith(DOT_PDB):
                pdb_paths.append(os.path.join(root, file))

    return pdb_paths


def schrod_prep(file_names):
    expected_files = set(map(lambda x: x + SCHROD_PREP_SUFF, file_names))

    shutil.copyfile(SCHROD_SCRIPT_PATH,
                    os.path.join(PDB_PREP_DIR, SCHROD_SCRIPT_PATH))

    os.chdir(PDB_PREP_DIR)
    subprocess.run(['bash', SCHROD_SCRIPT_PATH], stdout=subprocess.PIPE)

    res = {}

    while len(expected_files) > 0:
        print(expected_files)
        for path in os.listdir('.'):
            if path.endswith(SCHROD_PREP_SUFF) and path in expected_files:
                expected_files.remove(path)
                name = path[:-len(SCHROD_PREP_SUFF)]
                print(name)
                res[name] = os.path.join(PDB_PREP_DIR, path)
        time.sleep(1)

    os.chdir('..')

    return res


def prep_pdbs(last_epoch_name, epoch_name):
    pdb_paths = []

    for file in os.listdir(DB_PATH):
        dir_path = os.path.join(DB_PATH, file)
        if os.path.isdir(dir_path):
            pdb_paths += get_pdb_paths(dir_path, last_epoch_name)

    if not os.path.exists(PDB_PREP_DIR):
        os.mkdir(PDB_PREP_DIR)

    path_to_name = {}

    names = []

    for path in pdb_paths:
        name = path.replace('+', '_').replace('|', '_').split('/')
        name.reverse()
        name = '_'.join(name).replace(DOT_PDB, '') + DOT_PDB

        path_to_name[path] = name
        names.append(name)

        shutil.copyfile(path, os.path.join(PDB_PREP_DIR, name))

    path_to_prepped_schrod = schrod_prep(names)

    unprepped_paths = []

    for dir_path, name in path_to_name.items():
        if name not in path_to_prepped_schrod.keys():
            unprepped_paths.append(dir_path)
            continue
        path_to_prepped = path_to_prepped_schrod[name]

        new_path = dir_path.replace(last_epoch_name, epoch_name)

        if not os.path.exists(os.path.dirname(new_path)):
            os.makedirs(os.path.dirname(new_path))

        shutil.copyfile(path_to_prepped,
                        new_path)

    with open('unprepped_pdbs.csv', 'w') as f:
        for path in unprepped_paths:
            f.write(path + '\n')


if __name__ == '__main__':
    prep_pdbs(HETATMS_DELETED, PREPPED)
