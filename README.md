# abase

abase is a database of antibody-antigen complexes designed especially for
validation of algorithms of prediction of mutual conformation of
antibody and antigen.

Database consists of a **81** antibody-antigen complexs, for every one of which there
are both antibody and  antigen present in their *unbound* forms.

Information about complexes of the database and their unbound parts
can be found in the `abase_summary.csv`. Description of the summary 
file's columns can be found below.

Complexes themselves are stored in the `data` folder structure of
which is described below.

## abase_summary.csv

Summary file has following fields:

* **comp_name** — name of the given complex on terms of abase.
* **type** — type of complex. As of now there's only one type of 
complexes — **U:U**, what means that database contains both
unbound antibody and unbound antigen for complex. In the future there will
be two more types of complexes: **U:B**, **B:U**.