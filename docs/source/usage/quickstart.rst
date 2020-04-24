Quick Start
***************
If you installed PEPPA via PIP, run
::
  PEPPA --testunit

To download an examples/ folder which contains GFF files of four `Escherichia coli` ST131 genomes. 

If you installed PEPPA via git pull, the examples/ folder is already present in your PEPPA/ root folder. 

* Test PEPPA
::

  PEPPA -p examples/ST131 examples/*.gff.gz

* test PEPPA_parse
::

  PEPPA_parser -g examples/ST131.PEPPA.gff -s examples/PEPPA_out -m -t -c -a 95

