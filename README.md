## Installation ##
```
make FLORA

cd CSSR
make CSSR
cd ..

cd ViennaRNA
./configure --without-svm --without-swig --without-perl --without-python --without-doc-pdf --without-doc-html --without-doc --without-tutorial-pdf --without-tutorial --without-cla-pdf --without-cla --without-check --without-kinfold --without-forester --without-rnalocmin --without-rnaxplorer
make
```

## Usage ##

To see the help menu:
```
./FLORA.py
```

Run FLORA :
```
./FLORA.py input.fasta input.pdb output.pdb
```

To specify secondary structure in dot-bracket format:
```
./FLORA input.fasta input.pdb output.pdb input.dbn 3
```
Here, input.dbn should have only one line, where unpaired nucleotides are represented by '.', while paired nucleotides are represented by '(' or ')'.
