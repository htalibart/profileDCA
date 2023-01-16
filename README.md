# profileDCA

This package contains the following tools :

* *profileDCA_build* builds a profileDCA model from a sequence or a MSA (output is a folder which also contains positions matching between model and MSA / sequence files)
* *profileDCA_align* aligns profileDCA models and corresponding sequences / alignments
* *profileDCA_viz* allows you to visualize inferred profileDCA models
* *profileDCA_3Dviz* allows you to visualize the top N couplings of the inferred profileDCA models and whether the couplings are contacts in the 3D structure or not.

## Requirements

### Necessary requirements
profileDCA suite was developped with Python3.9 and requires the following packages, which will normally be automatically installed by setup.py:

* numpy (1.16.2)
* pandas (0.23.4)
* biopython (1.79)
* matplotlib (2.0.0)
* seaborn (0.9.0)

and, to build models from sequences, you will need HHsuite which you have to install and add to your path:

* HH-suite : https://github.com/soedinglab/hh-suite (tested with HHblits 3.0.3)
```
git clone https://github.com/soedinglab/hh-suite.git
mkdir -p hh-suite/build && cd hh-suite/build
cmake -DCMAKE_INSTALL_PREFIX=. ..
make -j 4 && make install
```
and export to your path :
```
export PATH="$(pwd)/bin:$(pwd)/scripts:$PATH"
```

### If you want to use profileDCA_3Dviz
To visualize predicted contacts, you also need to install:

* PyMOL : https://pymol.org/ (this package was developped for Pymol 2.1.0)

## Installation

### Compile the C++ solver library

```
cd profileDCA/profileDCA_align/apurva_cpp/
bash compile.bash
cd ..
```

### Install Python modules

```
pip3 install .
```


## Getting started
TODO


## License

This software is released under the terms of the GNU Affero General Public License v3.0 (GNU AGPLv3).

## Acknowledgments

The C++ ILP solver in profileDCA_align/apurva_cpp/ was designed by Wohlers, Andonov, Malod-Dognin, Klau and Yanev:
I. Wohlers, R. Andonov, G.W. Klau. DALIX: optimal DALI protein structure alignment. IEEE/ACM Trans Comput Biol Bioinform. 2013 Jan-Feb;10(1):26-36.
