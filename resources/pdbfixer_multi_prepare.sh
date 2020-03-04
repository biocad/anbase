#!/bin/bash
files=(*.pdb)
conda activate openmm
for f in $files
do
    echo "Prepping: $f" >> LOG.txt
	pdbfixer $f --replace-nonstandard --add-residues --output="$f.o.pdb"
	echo "Prepped: $f" >> LOG.txt
done