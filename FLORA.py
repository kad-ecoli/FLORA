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
os.environ['DFIRE_RNA_HOME']=os.path.join(rootdir,"dfire_rna")

def System(cmd):
    print(cmd)
    os.system(cmd)
    return

def check_consistent(inputFasta,inputPdb):
    sequence=''
    fp=open(inputFasta,'r')
    for line in fp.read().splitlines():
        if not line.startswith('>'):
            sequence+=line
    fp.close()
    sequence=sequence.lower()
    L=len(sequence)

    fp=open(inputPdb,'r')
    for line in fp.read().splitlines():
        if not line.startswith('ATOM  '):
            continue
        resSeq=int(line[22:26])
        nt    =line[19].lower()
        if resSeq<=0 or resSeq>L:
            print("ERROR! resSeq=%d out of range"%(resSeq))
            exit()
        elif sequence[resSeq-1]!=nt:
            print("ERROR! %s%d != %s"%(sequence[resSeq-1],resSeq,line))
    fp.close()
    return

if len(sys.argv)!=4:
    sys.stderr.write(docstring)
    exit()

inputFasta = sys.argv[1]
inputPdb   = sys.argv[2]
outputPdb  = sys.argv[3]

check_consistent(inputFasta,inputPdb)
if not os.path.isfile(outputPdb+".dbn"):
    System("%s/CSSR/CSSR %s %s %s.constraint"%(
        rootdir,inputFasta,inputPdb,outputPdb))
    System("%s/ViennaRNA/src/bin/RNAfold -i %s --noPS --constraint %s.constraint --enforceConstraint |head -2 |tail -1|grep -ohP '^\\S+'> %s.dbn"%(
        rootdir,inputFasta,outputPdb,outputPdb))
energy_list=[]
for m in [3,4,5]:
#for m in [0,1,2,3,4,5,6,7]:
    if not os.path.isfile("%s.%d"%(outputPdb,m)):
        System("%s/FLORA %s %s %s.%d %s.dbn %d"%(rootdir,inputFasta,inputPdb,outputPdb,m,outputPdb,m))
    if os.path.isfile("%s.%d"%(outputPdb,m)):
        cmd="%s/dfire_rna/bin/DFIRE_RNA %s.%d"%(rootdir,outputPdb,m)
        p=subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
        stdout,stderr=p.communicate()
        energy_list.append((float(stdout.decode().split()[1]),m))
energy_list.sort()
m=energy_list[0][1]
System("cp %s.%d %s"%(outputPdb,m,outputPdb))
