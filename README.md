# abase

abase is a database of antibody-antigen complexes designed especially for
validation of algorithms of prediction of mutual conformation of
antibody and antigen.

Database consists of a **86** antibody-antigen complexs, for every one of which there
are both antibody and  antigen present in their *unbound* forms.

Information about complexes of the database and their unbound parts
can be found in the `abase_summary.csv`. Description of the summary 
file's columns can be found below.

Complexes themselves are stored in the `data` folder structure of
which is described below.

## abase_summary.csv

Summary file has following fields:

* **comp_name** — name of the given complex in terms of abase.
* **type** — type of complex. As of now there's only one type of 
complexes — **U:U**, what means that database contains both
unbound antibody and unbound antigen for the 
complex. In the future there will
be two more types of complexes: **U:B** and **B:U**.
* **pdb_id_b** — PDB id of the *bound* complex.
* **ab_chain_ids_b** — ids of the 
antibody's chains in bound complex's structure.
* **ag_chain_ids_b** — ids of the 
antigen's chains in bound complex's structure.
* **ab_pdb_id_u** — PDB id of the unbound antigen.
* **ab_chain_ids_u** — ids of the 
unbound antibody's chains in its structure.
* **ag_pdb_id_u** — PDB id of the unbound antigen.
* **ag_chain_ids_u** — ids of the 
unbound antigen's chains in its structure.
* **ab_mismatches_cnt** — number of mismatches between unbound antibody's sequence and sequence of the antibody
from the bound complex.
* **ag_mismatches_cnt** — number of mismatches between unbound antigen's sequence and sequence of the antigen
from the bound complex.
* **small_molecules_message** — information about small molecules present in both
unbound antibody and unbound antigen. If equal to **NA**, 
all small molecules have no more than 7 heavy atoms. 
If all small molecules have no more than
15 heavy atoms, but some molecules have more than 7 heavy atoms, message 
**"small molecules with 7 < n_atoms <= 15 detected"** will be shown. If some
of the molecules have more than 15 heavy atoms, then message 
**"small molecules with n_atoms > 15 detected"** will be shown.
* **in_between_gaps_b** — number of gaps in the 
bound structures that are
fully contained in the radius of 15 Å of the interaction interface.
* **one_side_gaps_b** — number of gaps in the 
bound structures one end
of which is in the radius of 15 Å of the interaction interface and 
the other end is not.
* **long_gaps_b** — number of gaps in the bound structures
with length more than 15 amino acids.
* **total_gaps_b** — total number of gaps in the bound structures.
* **in_between_gaps_u** — number of gaps in the 
unbound structures that are
fully contained in the radius of 15 Å of the interaction interface.
* **one_side_gaps_u** — number of gaps in the 
unbound structures one end
of which is in the radius of 15 Å of the interaction interface and 
the other end is not.
* **long_gaps_u** — number of gaps in the unbound structures
with length more than 15 amino acids.
* **total_gaps_u** — total number of gaps in the unbound structures.
* **is_perfect** — if **True**, then **small_molecules_message** is equal to
**NA**, and both **in_between_gaps_u** and **one_side_gaps_u** are equal to
**0**. Otherwise is set to **False**.

## `data` folder

Every complex in the abase has its own subfolder in the `data`
folder. Complex's subfolder's name is equal to 
the corresponding **comp_name**.

Each **comp_name** folder contains following subfolders:

* `aligned` — folder that contains bound complex's structure and
the structures of its unbound components
structurally aligned to themselves in the bound complex.
* `hetatms_deleted` – folder that contains the same structures
as `aligned` folder with the HETATMs being removed
from the pdbs.
* `prepared_schrod` — folder that contains the same structures as
`hetatms_deleted` folder prepared with the use of Schrödinger's
*prepwizard*.
* `seqs` — folder that contains fasta-files with sequences of
the bound complex and its unbound structures.
* `annotation` — folder that contains annotation for the 
antibody in the bound complex and the antibody in its unbound
form.
* `constraints` — folder that contains constraints for the 
antigen in the bound complex and the antigen in its unbound
form.

Also, each **comp_name** folder contains file 
`alternative_candidates.csv`, that lists all the potential
candidates that can act as a pair of unbound antibody and 
unbound antigen for the given complex. 
This file has all the same fields
as `abase_summary.csv` with some changes:
* Field **candidate_name** is replaced with field **comp_name**.
* Field **is_perfect** is removed.
* Added fields **ca_rmsds**, **ncac_rmsds** and **all_atoms_rmsds** 
that contain corresponding rmsds from alternative candidate's 
interface to the main candidate's interface.

Subfolder of **comp_name** that is called `alternative_candidates`
contains candidates described in the `alternative_candidates.csv`.
Structure of the `alternative_candidates` folder is the same as
the structure of `data` with the exception that all its subfolders
are named after corresponding **candidate_name** fields. Also,
these subfolders don't contain info about alternative candidates,
because they represent these candidates themselves.

