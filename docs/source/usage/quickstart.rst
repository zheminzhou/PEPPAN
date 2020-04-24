Quick Start
***************
If you installed PEPPA via PIP, run
::

  PEPPA --testunit

To download an examples/ folder which contains GFF files of four *Escherichia coli* ST131 genomes. 

If you installed PEPPA via git pull, the examples/ folder is already present in your PEPPA/ root folder. 

* Test PEPPA

::

  PEPPA -p examples/ST131 -P examples/GCF_000010485.combined.gff.gz examples/*.gff.gz

* Test PEPPA_parse

::

  PEPPA_parser -g examples/ST131.PEPPA.gff -s examples/PEPPA_out -m -t -c -a 95

Key parameters for PEPPA
===========================
* GFFs

The main inputs for PEPPA. The nucleotide sequences can be integrated in the same file or in a separate file. Details see `Inputs <inputs.rst>`_

* -p or --prefix



* --min_cds
* --match_identity
* --pseudogene
