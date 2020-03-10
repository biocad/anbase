#!/bin/bash
files=(*$1.pdb)
maxtasks=50
ix=0
len=${#files[@]}
while [[ $ix -lt $len ]]
do
  echo "[$(date)] Started iteration $ix" >> LOG.txt
  tasks=$($SCHRODINGER/jobcontrol -list | wc -l)
  while [[ $tasks -le $maxtasks ]] && [[ $ix -lt $len ]]
  do
    f=${files[ix]}
    if test -f "$f.o.pdb"; then
        echo "[$(date)] Already prepared: $f" >> LOG.txt
    else
        $SCHRODINGER/utilities/prepwizard $(
            echo "-disulfides"                      # build the disulfide bridges
            echo "-fillloops"                       # fill missing residues
            echo "-fillsidechains"                  # fill side-chains where necessary
            echo "-mse"                             # convert selenomethionine residues to methionine
            echo "-noepik"                          # turn off epik, since it's used for small molecules and we don't have ones
            echo "-noimpref"                        # IMPORTANT: turn off minimization
            echo "-rehtreat"                        # IMPORTANT: delete existing hydrogen atoms and build new ones
            echo "-fasta_file $(realpath $f).fasta" # fasta file which is used to fill missing residues. IMPORTANT: has to be absolute
            echo "$f"
            echo "$f.o.pdb"
        )
        sleep 1
        tasks=$($SCHRODINGER/jobcontrol -list | wc -l)
        echo "[$(date)] Prepping: $f" >> LOG.txt
    fi
    let "ix += 1"
  done

  tasks=$($SCHRODINGER/jobcontrol -list | wc -l)
  while [[ $tasks -ne 1 ]]
  do
    sleep 5
    tasks=$($SCHRODINGER/jobcontrol -list | wc -l)
    echo "[$(date)] waiting" >> LOG.txt
  done
done
