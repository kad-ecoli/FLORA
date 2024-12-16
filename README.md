# FLORA: Full-atomic Loop Modeling of RNA #

Reconstruct missing nucleotides in the input RNA 3D structure.

## Installation ##

```
make FLORA

cd CSSR
make CSSR
cd ..

cd ViennaRNA
./configure --without-svm --without-swig --without-perl --without-python --without-doc-pdf --without-doc-html --without-doc --without-tutorial-pdf --without-tutorial --without-cla-pdf --without-cla --without-check --without-kinfold --without-forester --without-rnalocmin --without-rnaxplorer
make
cd ..

cd dfire_rna
make
cd ..
```

## Usage ##

### 1. Loop modeling ###

To see the help menu:
```
./FLORA.py
```

Run FLORA :
```
./FLORA.py input.fasta input.pdb output.pdb
```
Here, `input.fasta` is the FASTA format RNA sequence in lower case. `input.pdb` is the PDB format with missing nucleotides; the residue sequence number of `input.pdb` should be consistent with `input.fasta`. `input.pdb` can have missing atoms. `output.pdb` is the final full atomic model.

To specify secondary structure in dot-bracket format:
```
./FLORA input.fasta input.pdb output.pdb input.dbn 3
```
Here, input.dbn should have only one line, where unpaired nucleotides are represented by '.', while paired nucleotides are represented by '(' or ')'.

### 2. Full-atomic refinement ###

To just perform full-atomic structure refinement without loop modeling:
```bash
./Arena2 input.pdb output.pdb
```

## Citation ##
Junzhe Guoa, Lydia Freddolino and Chengxin Zhang (2025)
"FLORA: fast and accurate full-atomic loop modeling of RNA structures."
