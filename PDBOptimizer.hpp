#ifndef PDBOptimizer_HPP
#define PDBOptimizer_HPP 1

#include "PDBFiller.hpp"
#include "SSS.hpp"

using namespace std;

void write_config(char* filename_c, int option){
    //string filename_c=pdb_name+"_c.txt";
    ofstream fp(filename_c);
    string pdb_name=filename_c;
    pdb_name=pdb_name.substr(pdb_name.find_last_of('/')+1,pdb_name.find_last_of('_')-pdb_name.find_last_of('/')-1);
    fp<<"INPUTPDB   "<<"output_s/"<<pdb_name<<""<<".pdb"<<endl;
    fp<<"OUTPUTPDB  "<<"output_f/"<<pdb_name<<""<<".pdb"<<endl;
    fp<<"WRITEFREQ  "<<"100 "<<endl;
    fp<<"TRAJECTORY "<<"0       "<<endl;
    fp<<"VERBOSE    "<<"2       "<<endl;
    fp<<"NSTEPS     "<<"10000   "<<endl;
    fp<<"CUTOFF     "<<"12.0    "<<endl;
    fp<<"USEBORN    "<<"1       "<<endl;
    fp<<"ELECTR     "<<"1       "<<endl;
    fp<<"VDW        "<<"1       "<<endl;
    fp<<"HBONDS     "<<"1       "<<endl;
    fp<<"SSDETECT   "<<"1       "<<endl;
    fp<<"SSCONSTR   "<<"1       "<<endl;
    fp<<"NUMTHREADS  "<<"08      "<<endl;
    fp<<"BLOWGUARD  "<<"1       "<<endl;
    fp<<"POSRESTRAINTS "<<"1    "<<endl;
    fp<<"RESTRFILE "<<"restraints/"<<pdb_name<<"_r.txt"<<endl;
    fp.close();
}

void write_restraint(char* filename_dbn, char* filename_r, Model_F &fasta_in, ModelUnit &pdb_in, int option){
    Model_F dbn=read_dbn(filename_dbn,1);
    map<int,int> basepairs;
    char chain_ID=pdb_in.chains[0].chainID;
    for(int i=0;i<pdb_in.chains.size();i++){
        vector<int> pairif=pair_info(dbn.chains_F[i]);
        
        for(int i=0;i<pairif.size();i++){
            if(pairif[i]>0&&basepairs.find(pairif[i])==basepairs.end()) basepairs[i+1]=pairif[i];
        }
        for(auto iter:basepairs){
            cout<<iter.first<<"\t"<<iter.second<<endl;
        }
    }
    //string filename_r=pdb_name+"_r.txt";
    ofstream fp(filename_r);
    for(auto iter:basepairs){
        fp<<"BASEPAIR"<<"   "<<chain_ID<<"/"<<iter.first<<"   "<<chain_ID<<"/"<<iter.second<<"   "<<endl;
    }
    fp.close();
}

void LoopOptimize(ModelUnit &model_p,int chain,int gap_f,int gap_r){
}
#endif