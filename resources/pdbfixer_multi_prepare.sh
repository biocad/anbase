#!/bin/bash
files=(*.pdb)
for f in "${files[@]}"
do
    echo "Prepping: $f" >> LOG.txt
	$PDBFIXER $f --replace-nonstandard --add-residues --output="$f.o.pdb"
	echo "Prepped: $f" >> LOG.txt
done