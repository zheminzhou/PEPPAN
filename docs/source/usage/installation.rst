****************
Installation
****************
PEPPA has the following dependencies:

Required dependencies (You do not need to install any of these if you are based on a 64-bit Linux)

* `mmseqs2 <https://github.com/soedinglab/MMseqs2>`_
* `ncbi-blast+ <https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`_
* `diamond <https://github.com/bbuchfink/diamond>`_
* `rapidnj <https://birc.au.dk/software/rapidnj/>`_

Optional dependencies

* `fasttree2 <http://www.microbesonline.org/fasttree/#Install>`_ (For --orthology ml)

Installation of dependencies
----------------------------
* 64-bit Linux systems

Executables for all the dependencies are included as part of the PEPPA package. So by default, you do not need to install them. If the pre-compiled executables cannot run on your system, or if you want to use a different version of any dependency, please install the corresponding package in your system and point to the executables in your $PATH variable. 

The best way to install dependencies is via bio-conda. To do so, install `Anaconda <https://docs.anaconda.com/anaconda/install/>`_ or `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_. Then install bioconda and the dependencies.
::

  conda config --add channels defaults
  conda config --add channels conda-forge
  conda config --add channels bioconda
  conda install mmseqs2
  conda install blast
  conda install diamond
  conda install rapidnj
  conda install fasttree

You can check the presence of these dependencies with:
::

  command -v mmseqs blastn rapidnj diamond fasttree

Alternatively, you can always install each dependency by following the introductions in its link above. 

Install PEPPA
----------------------------
* Installing from PIP

::

  pip3 install bio-peppa


* Installing from source (advanced Linux users only)

Pull down the latest software from (https://github.com/zheminzhou/PEPPA).

Choose somewhere to put it, for example in your home directory (no root access required):
::

  cd $HOME
  git pull https://github.com/zheminzhou/PEPPA.git

And add the following lines to your $HOME/.bashrc file
::

  export PATH=$PATH:$HOME/PEPPA/


Test the installation
----------------------------

