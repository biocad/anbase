# anbase gathering pipeline

This directory contains all scripts that are necessary in order to
gather anbase from scratch or to update existing instance of the
database.

There're **9** stages of the anbase's gathering, each one of which is accompanied
by a corresponding script:

1. `fetch_unbound_data.py` — script that given relevant `sabdab_summary_all.tsv`,
that can be downloaded from [here](http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/search/?all=true#downloads), finds candidates for unbound parts
of SAbDab's complexes.

2. `filter_unbound_data.py` — script that determines whether unbound candidates
for structures that were found as the result of the previous script's
work are indeed unbound.

3. `process_unbound_data.py` — script that downloads unbound structures for
the complexes, structurally aligns these structures on their corresponding
parts in the complexes and stores all structures and complexes in the
`data` folder. Script also forms .csv-table `db_info.csv`, that
contains information about complexes from SAbDab for which there
were found unbound parts.

4. `prepper.py` — script that fixes broken structures and restores missing hydrogens.
Since all structures contained in anbase are downloaded from PDB, the chances
are they have serious flaws in them: missing atoms, missing residues, etc.
In order to fix flawed structures a particular instrument is needed.
`prepper.py` uses either a 
[Schrödinger PrepWizard](https://www.schrodinger.com/protein-preparation-wizard) 
or a [PDBFixer](https://github.com/openmm/pdbfixer) as
such instrument. If you don't have any of them installed in your system,
script won't work.

5. `gap_stats_collector.py` — script that checks whether there are any gaps — absent residues —
in the anbase's structures after the work of `prepper.py`.

6. `duplicates_analyser.py` — script that performs elimination of duplicates
in the database. Important requirement to run this script is the presence
of annotation for the antibodies in the database. One of the ways to obtain 
annotation for the anbase's antibodies is to use 
[this script](https://github.com/biocad/hedge/blob/alexkanerus/develop/app/AnnotateAbase.hs) (link leads to
*private* BIOCAD's repository).

7. `generate_constraints.py` — script that generates docking constraints
for the given complex. During anbase's gathering it's used to generate
constraints for the native pose of the epitope.

8. `form_anbase_summary.py` — script that finalizes all results of the 
previous scripts:

    * it creates `anbase_summary.csv` in the format described [here](https://github.com/biocad/anbase/blob/master/README.md);
    
    * it creates directory that is structured like [`data`-directory](https://github.com/biocad/anbase/tree/master/data)
    from the root of the anbase repository.

9. `aligner.py` — script that aligns alternative candidates on the complexes
from the main ones. Works on the directory produced by the `form_anbase_summary.py`
script.