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
    $SCHRODINGER/utilities/prepwizard -disulfides -fillloops -fillsidechains -mse -noepik -noimpref -fasta_file "$f.fasta" -rehtreat "$f" "$f.o.pdb"
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