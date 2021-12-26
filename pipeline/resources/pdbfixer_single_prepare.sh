#!/bin/bash
input_file=$1

if test -f "$input_file.o.pdb"; then
    echo "Already prepared: $input_file" >> LOG.txt
else
    echo "[$(date)] Prepping: $input_file" >> LOG.txt
    $PDBFIXER $input_file --replace-nonstandard --add-residues --output="$input_file.o.pdb"
    echo "[$(date)] Prepped: $input_file" >> LOG.txt
fi