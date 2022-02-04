# anbase

anbase is a database of antibody-antigen complexes designed especially for
validation of algorithms of prediction of mutual conformation of
antibody and antigen.

Database consists of a **570** antibody-antigen complexs: **75** of which are *unbound-antibody*:*unbound-antigen*, and **495** are *bound-antibody:unbound-antigen*.

Information about complexes of the database and their unbound parts
can be found in the `anbase_summary.csv`. Description of the summary 
file's columns can be found below.

Complexes themselves are stored in the `data` folder structure of
which is described below.

## anbase_summary.csv

Summary file has following fields:

* **comp_name** — name of the given complex in terms of anbase.
* **type** — type of complex. There are several types of complexes. **U:U** type means that database contains both
unbound antibody and unbound antigen for the 
complex. **B:U** — bound antibody and unbound antigen.
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

Every complex in the anbase has its own subfolder in the `data`
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
* `epitope` — folder that contains constraints for the 
antigen's epitope in the bound complex and the antigen's epitope 
in its unbound form.

Each **comp_name** may contain several alternative unbound versions. All versions are in folder named with numbers starting with 0.
