#ifndef ADJUSTPOSITION_HPP
#define ADJUSTPOSITION_HPP

#include "cssr.hpp"
#include "BondLengths.hpp"
#include "BaseConformation.hpp"
#include "AtomicClashes.hpp"
#include "Chirality.hpp"

void adjust_position(ModelUnit &pdb_in,vector<vector<size_t> > &res_str_vec,
                     vector<pair<double,vector<size_t> > > &bp_vec,
                     map<string, map<string,vector<double> > > &ideal_rna,
                     const int typing=1,const int option=0,double tolerance=0.1){
    vector<vector<int> >moveable_mat;
    check_moveable(pdb_in, moveable_mat);
    for (int c=0; c<pdb_in.chains.size(); c++){
        int length = pdb_in.chains[c].residues.size();
        //int iterations = round(0.08*length+500);
        int iterations = 1000;
        for (int t=0; t<iterations; t++){
            //cout<<"A"<<endl;
    	    size_t moved = fix_bond_lengths(pdb_in, tolerance);  
            moved+=fix_chirality(pdb_in);
            //cout<<"B"<<endl;
            if(option>1){
		//moved+=fix_chirality(pdb_in);
                if (t%10==0)
                {
                    size_t bp;
                    for (bp=0;bp<res_str_vec.size();bp++) res_str_vec[bp].clear();
                    for (bp=0;bp<res_str_vec.size();bp++) res_str_vec[bp].shrink_to_fit();
                    res_str_vec.clear();
                    res_str_vec.shrink_to_fit();
                    for (bp=0;bp<bp_vec.size();bp++) bp_vec[bp].second.clear();
                    for (bp=0;bp<bp_vec.size();bp++) bp_vec[bp].second.shrink_to_fit();
                    bp_vec.clear();
                    bp_vec.shrink_to_fit();
                    cssr(pdb_in, res_str_vec, bp_vec);
                    filter_bp(bp_vec);
                }
                moved+=fixBaseConformation(pdb_in, ideal_rna, bp_vec);
            }
            if (t+1<iterations) moved+=fix_clashes(pdb_in, moveable_mat);
            else moved+=fix_bond_lengths(pdb_in, tolerance);  
            //cout<<"C"<<endl;
                // if no atoms are moved, immediately break out of this for loop
        	if (moved==0) break;
            if(typing) cout<<"t="<<t<<" moved="<<moved<<endl;
        }
    }
    vector<vector<size_t> >().swap(res_str_vec);
    vector<pair<double,vector<size_t> > > ().swap(bp_vec);
}

#endif
