#!/bin/bash

input_file=$1

echo "schrod_single_prepare.sh. input_file: $input_file"

$SCHRODINGER/utilities/prepwizard $(
    echo "-disulfides"                      # build the disulfide bridges
    echo "-fillloops"                       # fill missing residues
    echo "-fillsidechains"                  # fill side-chains where necessary
    echo "-mse"                             # convert selenomethionine residues to methionine
    echo "-noepik"                          # turn off epik, since it's used for small molecules and we don't have ones
    echo "-noimpref"                        # IMPORTANT: turn off minimization
    echo "-rehtreat"                        # IMPORTANT: delete existing hydrogen atoms and build new ones
    echo "-fasta_file $(realpath $input_file).fasta" # fasta file which is used to fill missing residues. IMPORTANT: has to be absolute
    echo "-WAIT"                            # do it synchronously
    echo "$input_file"
    echo "$input_file.o.pdb"
)