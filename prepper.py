import multiprocessing
import os
import shutil
import subprocess
import time

DB_PATH = 'data'
DOT_PDB = '.pdb'

HETATMS_DELETED = 'hetatms_deleted'

PDB_PREP_DIR = 'pdb_prep'
SCHROD_SCRIPT_PATH = 'resources/schrod_multi_prepare.sh'

PREP_SUFF = '.o.pdb'

PREPPED = 'prepared'

SCHROD = 'schrod'
PDB_FIXER = 'PDBFixer'

SEQS = 'seqs'

DOT_FASTA = '.fasta'


def get_pdb_paths(dir_path, prev_epoch):
    pdb_paths = []

    for root, _, files in os.walk(dir_path):
        for file in files:
            if prev_epoch in root and file.endswith(DOT_PDB):
                pdb_paths.append(os.path.join(root, file))

    return pdb_paths


def await_expected_files(expected_files, tmp_dir):
    res = {}

    waiting_time = 10 * 60

    while len(expected_files) > 0:
        print(expected_files)
        for path in os.listdir('.'):
            if path.endswith(PREP_SUFF) and path in expected_files:
                expected_files.remove(path)
                name = path[:-len(PREP_SUFF)]
                print(name)
                res[name] = os.path.join(tmp_dir, path)
        time.sleep(1)
        waiting_time -= 1

        if waiting_time == 0:
            break

    return res


def schrod_prep(file_names, tmp_dir):
    expected_files = set(map(lambda x: x + PREP_SUFF, file_names))

    script_name = os.path.basename(SCHROD_SCRIPT_PATH)

    shutil.copyfile(SCHROD_SCRIPT_PATH,
                    os.path.join(tmp_dir, script_name))

    os.chdir(tmp_dir)
    subprocess.run(['bash', script_name], stdout=subprocess.PIPE)

    res = await_expected_files(expected_files, tmp_dir)

    os.chdir('..')

    return res


def run_pdb_fixer(name):
    print('Prepping:', name, flush=True)
    subprocess.run(
        ['pdbfixer', name, '--replace-nonstandard', '--add-residues',
         '--output={}{}'.format(name, PREP_SUFF)], stdout=subprocess.PIPE)
    print('Prepped:', name, flush=True)


def pdb_fixer_prep(file_names, tmp_dir):
    expected_files = set(map(lambda x: x + PREP_SUFF, file_names))

    os.chdir(tmp_dir)

    pool = multiprocessing.Pool()
    r = pool.map_async(run_pdb_fixer, file_names)
    r.wait()

    res = await_expected_files(expected_files, tmp_dir)

    os.chdir('..')

    return res


def prep_pdbs(last_epoch_name, epoch_name, db_path, mode, tmp_dir):
    pdbs_to_copy = []

    for file in os.listdir(db_path):
        dir_path = os.path.join(db_path, file)
        if os.path.isdir(dir_path):
            pdbs_to_copy += get_pdb_paths(dir_path, last_epoch_name)

    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)

    path_to_name = {}

    names = []

    for path in pdbs_to_copy:
        name = path.replace('+', '_').replace('|', '_').split('/')
        name.reverse()
        name = '_'.join(name).replace(DOT_PDB, '') + DOT_PDB

        path_to_name[path] = name
        names.append(name)

        shutil.copyfile(path.replace(last_epoch_name, SEQS).
                        replace(DOT_PDB, DOT_FASTA),
                        os.path.join(tmp_dir, name) + DOT_FASTA)
        shutil.copyfile(path, os.path.join(tmp_dir, name))

    path_to_prepped_schrod = schrod_prep(names, tmp_dir) if mode == SCHROD \
        else pdb_fixer_prep(names, tmp_dir)

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
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('--db', default=DB_PATH, dest='db', metavar='DB',
                      help='Path to database [default: {}]'.format(DB_PATH))
    parser.add_option('--mode', default=SCHROD, dest='mode', metavar='MODE',
                      help='{} for preparation using schrodinger. '
                           '{} for preparation using PDBFixer. [default: {}]'.
                      format(SCHROD, PDB_FIXER, SCHROD))
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
    options, _ = parser.parse_args()

    prep_pdbs(options.prev_epoch, options.cur_epoch, options.db, options.mode,
              options.tmp_dir)
