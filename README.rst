Protein DataBank file parsing utilities
=======================================

This project provides a small set of utilities to parse PDB files. It is not
written using `bioptyhon <http://biopython.org/>`_ and thus selfcontained.

Reading PDB-files into memory, selecting chains, renumbering/removing atoms and
extracting distance matrices are the currently supported operations. See
examples below.

Usage examples
--------------

Extracting a distance matrix, consider only CA atoms;

   pdbutil_remove_atoms -keep CA <pdbfile> | pdbutil_distances

Test installation
-----------------

Install this package under a local user::

   # Install editable under local user
   pip install --user -e .;
   # Export ~/.local into PATH, if not there already
   export PATH=${HOME}/.local/bin:${PATH};