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
GFFs
--------------------

The main inputs for PEPPA. The nucleotide sequences can be integrated in the same file or in a separate file. Details see `Inputs <inputs.rst>`_

-p or --prefix
--------------------

The prefix for the output files, as well as intermediate files. The final outputs are:

1. <prefix>.PEPPA.gff
2. <prefix>.allele.fasta

Details see `Outputs <outputs.rst>`_.

--match_identity
--------------------
match_identity controls the minimal identity of an alignment to be considered in the pan-genome construction. PEPPA is insensitive to the value of this parameter, as long as it is low enough to cover the genetic diversity of the genomes. However, reducing this value will increase the run time of the program. By default, PEPPA accepts a minimum match identity of 65%. Use a lower value if you run PEPPA between genomes from different genera.

--min_cds
--------------------
This controls the minimal size of a gene to be considered in the final pan-genome. Short genes are less reliable and also hard to detect in the similarity-based search. By default, PEPPA ignores genes with <= 150 bps in their sizes. However, you can always reduce the minimum size of the accepted genes if trust the original annotations.

--pseudogene
--------------------
A coding region in the genome is assigned as a pseudogene, if its size is shorter than other orthologous genes by a certain proportion. By default, PEPPA sets --pseudogene 0.8, therefore any gene that is >=80% of the representative gene will be assign as an "intact" gene, otherwise a "pseudogene". 

See `Parameters <parameters.rst>`_ for description of other paramters implemented in PEPPA. 
