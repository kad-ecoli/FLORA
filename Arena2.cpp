#include <iostream>
#include <iomanip>
#include <fstream>
#include "MissingRNAatom.hpp"
#include "cssr.hpp"
#include "AdjustPosition.hpp"
#include "PDBParser.hpp"
#include "SSS.hpp"
//#include "PDBOptimizer.hpp"


const char* helptxt=""
"Arena2 input.pdb output.pdb\n"
"    Find and fill missing residues and atoms in input.pdb\n"
"    Output filled model to output.pdb\n"
;


int main(int argc,char **argv){
    /*input*/
    if (argc<=2)
    {
        cerr<<helptxt;
        return 0;
    }
    string infile_PDB=argv[1];
    string outfile   =argv[2];
    
    /*build data structures*/
    ModelUnit pdb_in=read_pdb_structure(infile_PDB.c_str(),2,0);

    /*CSSR and Fill missing Atoms*/
    vector<vector<size_t> >res_str_vec;
    vector<pair<double,vector<size_t> > > bp_vec;
    cssr(pdb_in, res_str_vec, bp_vec);
    filter_bp(bp_vec);
    map<string, map<string,vector<double> > > ideal_rna=parse_ideal_rna();
    MissingRNAatom(pdb_in,ideal_rna,bp_vec,5);
    standardize_pdb_order(pdb_in);
    adjust_position(pdb_in,res_str_vec,bp_vec,ideal_rna);
    size_t a,r,c;
    for (c=0;c<pdb_in.chains.size();c++)
        for (r=0;r<pdb_in.chains[c].residues.size();r++)
            for (a=0;a<pdb_in.chains[c].residues[r].atoms.size();a++)
                if (!pdb_in.chains[c].residues[r].atoms[a].movable)
                    pdb_in.chains[c].residues[r].atoms[a].movable=1;
    adjust_position(pdb_in,res_str_vec,bp_vec,ideal_rna);
    
    /*output and clear the maps*/
    map<string, map<string,vector<double> > >().swap(ideal_rna);
    modify_RNA_structure(pdb_in);
    write_pdb_structure(outfile.c_str(),pdb_in);
    vector<ChainUnit>().swap(pdb_in.chains);
    return 0;
}
