#!/bin/bash
files=(*$1.pdb)
for f in "${files[@]}"
do
    if test -f "$f.o.pdb"; then
        echo "Already prepared: $f" >> LOG.txt
    else
        echo "Prepping: $f" >> LOG.txt
        $PDBFIXER $f --replace-nonstandard --add-residues --output="$f.o.pdb"
        echo "Prepped: $f" >> LOG.txt
    fi
done