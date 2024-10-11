#include <iostream>
#include <iomanip>
#include <fstream>
#include "MissingRNAatom.hpp"
#include "cssr.hpp"
#include "AdjustPosition.hpp"
#include "PDBParser.hpp"
#include "PDBFiller.hpp"
#include "SSS.hpp"
//#include "PDBOptimizer.hpp"


const char* helptxt=""
"FLORA sequence.fasta input.pdb output.pdb pair.dbn [option] [tolerance]\n"
"    Find and fill missing residues and atoms in input.pdb\n"
"    Output filled model to output.pdb\n"
"option:\n"
;


int main(int argc,char **argv){
    /*input*/
    if (argc<4){
        cerr<<helptxt;
        return 0;
    }
    string infile_fasta=argv[1];
    string infile_PDB=argv[2];
    //cout<<infile_PDB<<endl;
    string outfile=argv[3];
    string dbn=argv[4];
    int option=(argc>5)?atoi(argv[5]):3;
    double tolerance =(argc>6)?atof(argv[6]):0.1; // default tolerance is 10%
    
    /*build data structures*/
    Model_F fasta_in=read_fasta_structure(argv[1]);
    ModelUnit pdb_in=read_pdb_structure(argv[2],2,0);

    /*examine infiles*/

    if(compare_fasta_pdb(fasta_in,pdb_in)) return 0;
    
    /*CSSR and Fill missing Atoms*/
    
    fill_residue(fasta_in,pdb_in);//fill residue name for cssr() and MissingRNAatom()
    vector<vector<size_t> >res_str_vec;
    vector<pair<double,vector<size_t> > > bp_vec;
    cssr(pdb_in, res_str_vec, bp_vec);
    filter_bp(bp_vec);
    map<string, map<string,vector<double> > > ideal_rna=parse_ideal_rna();
    MissingRNAatom(pdb_in,ideal_rna,bp_vec,5);
    standardize_pdb_order(pdb_in);
    adjust_position(pdb_in,res_str_vec,bp_vec,ideal_rna);
    
    vector<sizet_sizet_vc> gapif_s=put_single_pairs(argv[4],fasta_in,pdb_in,ideal_rna);
    MissingRNAatom(pdb_in,ideal_rna,bp_vec,5);
    standardize_pdb_order(pdb_in);
    adjust_position(pdb_in,res_str_vec,bp_vec,ideal_rna,1);

    /*fill the gap*/
   
    choose_fill_strategy(fasta_in,pdb_in,ideal_rna,bp_vec,res_str_vec,option);
    
    put_double_pairs(argv[4],fasta_in,pdb_in,gapif_s,bp_vec,ideal_rna);
    MissingRNAatom(pdb_in,ideal_rna,bp_vec,5);
    standardize_pdb_order(pdb_in);
    adjust_position(pdb_in,res_str_vec,bp_vec,ideal_rna,1);
    //write_pdb_structure(outfile.c_str(),pdb_in);

    choose_fill_strategy(fasta_in,pdb_in,ideal_rna,bp_vec,res_str_vec,option);
    
    //write_pdb_structure(outfile.c_str(),pdb_in);

    /*CSSR again*/
    cssr(pdb_in, res_str_vec, bp_vec);
    filter_bp(bp_vec);
    MissingRNAatom(pdb_in,ideal_rna,bp_vec,5);
    //write_pdb_structure(outfile.c_str(),pdb_in);
    standardize_pdb_order(pdb_in);

    adjust_position(pdb_in,res_str_vec,bp_vec,ideal_rna,1,2);
    //write_pdb_structure(outfile.c_str(),pdb_in);
    /*output and clear the maps*/
    map<string, map<string,vector<double> > >().swap(ideal_rna);
    vector<sizet_sizet_vc>().swap(gapif_s);
    modify_RNA_structure(pdb_in);
    write_pdb_structure(outfile.c_str(),pdb_in);
    vector<ChainUnit>().swap(pdb_in.chains);
    
    
    return 0;
}
