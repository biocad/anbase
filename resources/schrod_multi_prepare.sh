#!/bin/bash
files=(*.pdb)
ix=0
iy=0
len=${#files[@]}
while [[ $ix -le $len ]]
do
  let 'iy = ix + 50'
  ctime=$(date)
  echo "[$ctime] Started iteration $ix" >> LOG.txt
  curbatch=${files[@]:$ix:50}
  for f in  $curbatch
  do
    $SCHRODINGER/utilities/prepwizard $(
        echo "-disulfides"          # build the disulfide bridges
        echo "-fillloops"           # fill missing residues
        echo "-fillsidechains"      # fill side-chains where necessary
        echo "-mse"                 # convert selenomethionine residues to methionine
        echo "-noepik"              # turn off epik, since it's used for small molecules and we don't have ones
        echo "-noimpref"            # IMPORTANT: turn off minimization
        echo "-rehtreat"            # IMPORTANT: delete existing hydrogen atoms and build new ones
        echo "-fasta_file $(realpath $f).fasta" # fasta file which is used to fill missing residues. IMPORTANT: has to be absolute
        echo "$f"
        echo "$f.o.pdb"
    )
  done
  let 'ix += 50'
  sleep 5
  tasks=$($SCHRODINGER/jobcontrol -list | wc -l)
  while [[ $tasks -gt "1" ]]
  do
    sleep 5
    ctime=$(date)
    echo "[$ctime] Waiting for $tasks tasks to terminate. Iteration $ix" >> LOG.txt
    tasks=$($SCHRODINGER/jobcontrol -list | wc -l)
  done
done