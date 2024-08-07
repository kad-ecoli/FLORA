#!/usr/bin/python3
docstring='''
FLORA.py input.fasta input.pdb output.pdb

Input:
    input.fasta - fasta format sequence
    input.pdb   - input pdb format with partial structure

Output:
    output.pdb  - full atomic final model
'''
import sys
import os
import subprocess
from math import pi, sqrt, acos

rootdir=os.path.dirname(os.path.abspath(__file__))

def System(cmd):
    print(cmd)
    os.system(cmd)
    return

if len(sys.argv)!=4:
    sys.stderr.write(docstring)
    exit()

inputFasta = sys.argv[1]
inputPdb   = sys.argv[2]
outputPdb  = sys.argv[3]

System("%s/CSSR/CSSR %s %s %s.constraint"%(
        rootdir,inputFasta,inputPdb,outputPdb))
System("%s/ViennaRNA/src/bin/RNAfold -i %s --noPS --constraint %s.constraint --enforceConstraint |head -2 |tail -1|grep -ohP '^\\S+'> %s.dbn"%(
        rootdir,inputFasta,outputPdb,outputPdb))
System("%s/FLORA %s %s %s %s.dbn 3"%(rootdir,inputFasta,inputPdb,outputPdb,outputPdb))
