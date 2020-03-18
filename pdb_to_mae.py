import os
import subprocess

HETATMS_DELETED = 'hetatms_deleted'
DB_PATH = 'data'

DOT_PDB = '.pdb'
DOT_MAE = '.mae'


def convert_pdb_paths(dir_path, cur_epoch):
    for root, _, files in os.walk(dir_path):
        for file in files:
            if cur_epoch in root and file.endswith(DOT_PDB):
                path_to_pdb = os.path.join(os.path.abspath(root), file)
                path_to_mae = path_to_pdb[:-4] + DOT_MAE

                if os.path.exists(path_to_mae):
                    continue
                    
                command = '$SCHRODINGER/utilities/structconvert ' \
                          '-ipdb \'{}\' -omae \'{}\''.format(path_to_pdb, path_to_mae)
                subprocess.call(command, stdout=subprocess.PIPE, shell=True)


def convert_pdbs(db_path, cur_epoch):
    for file in os.listdir(db_path):
        dir_path = os.path.join(db_path, file)
        if os.path.isdir(dir_path):
            print(file, flush=True)
            convert_pdb_paths(dir_path, cur_epoch)


if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('--db', default=DB_PATH, dest='db', metavar='DB',
                      help='Path to database [default: {}]'.format(DB_PATH))
    parser.add_option('--cur-epoch', default=HETATMS_DELETED,
                      dest='cur_epoch', metavar='CUR_EPOCH',
                      help='Name of the epoch structures from which will be '
                           'converted to mae. [default: {}]'.format(
                          HETATMS_DELETED))
    options, _ = parser.parse_args()

    convert_pdbs(options.db, options.cur_epoch)
