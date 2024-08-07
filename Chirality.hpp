#ifndef Chirality_HPP
#define Chirality_HPP 1

#include <math.h>
#include <stdbool.h>
#include <vector>
#include <cstdlib>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <map>

#include "GeometryTools.hpp"

using namespace std;

//224.rna.pdb A A   5

double P_chirality[4][3]={
    {4.598, -8.226, -15.542},//OP1
    {4.155, -6.085, -14.310},//OP2
    {4.074, -7.563, -14.327},//P
    {4.807, -8.157, -13.034}//O3' n-1
};

double C2_chirality[4][3]={
    {-0.707, -8.471, -12.069},//C1'
    {-0.157, -8.816, -14.332},//C3'
    {-1.368, -8.816, -13.400},//C2'
    {-2.138, -9.993, -13.260}//O2'
};

double C4_chirality[4][3]={
    {-0.157, -8.816, -14.332},//C3'
    {2.293, -9.405, -13.779},//C5'
    {0.579, -9.069, -12.111},//O4'
    {0.825, -9.588, -13.450}//C4'
};

double C3_chirality[4][3]={
    {-1.368, -8.816, -13.400},//C2'
    {0.825, -9.588, -13.450},//C4'
    {-0.361, -9.460, -15.582},//O3'
    {-0.157, -8.816, -14.332}//C3'
};

//224.rna.pdb G A   6

double C1_N9_chirality[4][3]={
    {-4.412, -7.943, -14.659},//O4'
    {-5.913, -6.681, -15.948},//C2'
    {-4.252, -5.596, -14.387},//N9
    {-5.171, -6.746, -14.617}//C1'
};

//224.rna.pdb C A   4

double C1_N1_chirality[4][3]={
    {5.386, -7.318, -9.563},//O4'
    {3.612, -8.157, -10.852},//C2'
    {3.318, -6.195, -9.291},//N1
    {3.981, -7.510, -9.521}//C1'
};

inline int find_atom_in_residue(ResidueUnit &residue,string atom_name){
    for(int i=0;i<residue.atoms.size();i++){
        if(residue.atoms[i].name==atom_name) return i;
    }
    return -1;
}

bool OnSameLine(vector<vector<double> > &atoms){
    
    if(atoms.size()>3) return false;
    if(atoms.size()<3) return true;
    for(int i=0;i<3;i++){
        if(atoms[i].size()<3) return true;
    }
    
    vector<double> delta1(3,0);
    vector<double> delta2(3,0);
    for(int i=0;i<3;i++) delta1[i]=atoms[0][i]-atoms[1][i];
    for(int i=0;i<3;i++) delta2[i]=atoms[1][i]-atoms[2][i];
    
    double a2,b2,ab,alpha,dev;
    a2=delta1[0]*delta1[0]+delta1[1]*delta1[1]+delta1[2]*delta1[2];
    b2=delta2[0]*delta2[0]+delta2[1]*delta2[1]+delta2[2]*delta2[2];
    ab=delta1[0]*delta2[0]+delta1[1]*delta2[1]+delta1[2]*delta2[2];
    
    vector<double>().swap(delta1);
    vector<double>().swap(delta2);
    
    if(sqrt(a2*b2)<Extra) return true;
    alpha=ab/sqrt(a2*b2);
    dev=abs(1-abs(alpha));
    if(dev<Extra) return true;//const double Extra=1.0e-4;
    else return false;
}

int check_carbon_chirality(ResidueUnit &residue,string atom_name){
    //exclude the cases where atoms cannot be moved
    if((atom_name==" C2'"&&!residue.atoms[find_atom_in_residue(residue," O2'")].movable)||
       (atom_name==" C4'"&&!residue.atoms[find_atom_in_residue(residue," C4'")].movable)||
       (atom_name==" C3'"&&!residue.atoms[find_atom_in_residue(residue," C3'")].movable)||
       (atom_name==" C1'"&&!residue.atoms[find_atom_in_residue(residue," C1'")].movable)) return 0;
    vector<double> tmp(3,0);
    vector<vector<double> > ideal_true_list;
    vector<vector<double> > real_list;
    ideal_true_list.assign(4,tmp);
    real_list.assign(4,tmp);
    if(atom_name==" C2'"){
        for(int j=0;j<4;j++){
            for(int k=0;k<3;k++){
                ideal_true_list[j][k]=C2_chirality[j][k];
            }
        }
        
        for(int j=0;j<3;j++) real_list[0][j]=residue.atoms[find_atom_in_residue(residue," C1'")].xyz[j];
        for(int j=0;j<3;j++) real_list[1][j]=residue.atoms[find_atom_in_residue(residue," C3'")].xyz[j];
        for(int j=0;j<3;j++) real_list[2][j]=residue.atoms[find_atom_in_residue(residue," C2'")].xyz[j];
        for(int j=0;j<3;j++) real_list[3][j]=residue.atoms[find_atom_in_residue(residue," O2'")].xyz[j];

    }else if(atom_name==" C4'"){
        for(int j=0;j<4;j++){
            for(int k=0;k<3;k++){
                ideal_true_list[j][k]=C4_chirality[j][k];
            }
        }

        for(int j=0;j<3;j++) real_list[0][j]=residue.atoms[find_atom_in_residue(residue," C3'")].xyz[j];
        for(int j=0;j<3;j++) real_list[1][j]=residue.atoms[find_atom_in_residue(residue," C5'")].xyz[j];
        for(int j=0;j<3;j++) real_list[2][j]=residue.atoms[find_atom_in_residue(residue," O4'")].xyz[j];
        for(int j=0;j<3;j++) real_list[3][j]=residue.atoms[find_atom_in_residue(residue," C4'")].xyz[j];

    }else if(atom_name==" C3'"){
        for(int j=0;j<4;j++){
            for(int k=0;k<3;k++){
                ideal_true_list[j][k]=C3_chirality[j][k];
            }
        }
        for(int j=0;j<3;j++) real_list[0][j]=residue.atoms[find_atom_in_residue(residue," C2'")].xyz[j];
        for(int j=0;j<3;j++) real_list[1][j]=residue.atoms[find_atom_in_residue(residue," C4'")].xyz[j];
        for(int j=0;j<3;j++) real_list[2][j]=residue.atoms[find_atom_in_residue(residue," O3'")].xyz[j];
        for(int j=0;j<3;j++) real_list[3][j]=residue.atoms[find_atom_in_residue(residue," C3'")].xyz[j];

    }else if(atom_name==" C1'"&&(residue.resn=="  A"||residue.resn=="  G")){
        for(int j=0;j<4;j++){
            for(int k=0;k<3;k++){
                ideal_true_list[j][k]=C1_N9_chirality[j][k];
            }
        }

        for(int j=0;j<3;j++) real_list[0][j]=residue.atoms[find_atom_in_residue(residue," O4'")].xyz[j];
        for(int j=0;j<3;j++) real_list[1][j]=residue.atoms[find_atom_in_residue(residue," C2'")].xyz[j];
        for(int j=0;j<3;j++) real_list[2][j]=residue.atoms[find_atom_in_residue(residue," N9 ")].xyz[j];
        for(int j=0;j<3;j++) real_list[3][j]=residue.atoms[find_atom_in_residue(residue," C1'")].xyz[j];

    }else if(atom_name==" C1'"&&(residue.resn=="  U"||residue.resn=="  C")){
        for(int j=0;j<4;j++){
            for(int k=0;k<3;k++){
                ideal_true_list[j][k]=C1_N1_chirality[j][k];
            }
        }

        for(int j=0;j<3;j++) real_list[0][j]=residue.atoms[find_atom_in_residue(residue," O4'")].xyz[j];
        for(int j=0;j<3;j++) real_list[1][j]=residue.atoms[find_atom_in_residue(residue," C2'")].xyz[j];
        for(int j=0;j<3;j++) real_list[2][j]=residue.atoms[find_atom_in_residue(residue," N1 ")].xyz[j];
        for(int j=0;j<3;j++) real_list[3][j]=residue.atoms[find_atom_in_residue(residue," C1'")].xyz[j];

    }
    if(ideal_true_list.size()==0||real_list.size()==0){
        vector<vector<double> >().swap(ideal_true_list);
        vector<vector<double> >().swap(real_list);
        return 0;
    }
    
    vector<vector<double> > RotMatix;
    vector<double> TranVect;

    RotateCoor(ideal_true_list,real_list,RotMatix,TranVect);
    double SD_true=0;
    for(int i=0;i<4;i++){
        ChangeCoor(ideal_true_list[i],RotMatix,TranVect,tmp);
        SD_true+=Points2Distance2(real_list[i],tmp);
    }
    vector<vector<double> >().swap(RotMatix);
    vector<double>().swap(TranVect);

    ideal_true_list[0].swap(ideal_true_list[1]);
    RotateCoor(ideal_true_list,real_list,RotMatix,TranVect);
    double SD_false=0;
    for(int i=0;i<4;i++){
        ChangeCoor(ideal_true_list[i],RotMatix,TranVect,tmp);
        SD_false+=Points2Distance2(real_list[i],tmp);
    }
    vector<vector<double> >().swap(RotMatix);
    vector<double>().swap(TranVect);
    if(SD_true<SD_false){
        vector<vector<double> >().swap(ideal_true_list);
        vector<vector<double> >().swap(real_list);
        return 0;
    }
    tmp.swap(ideal_true_list[3]);
    ideal_true_list.pop_back();
    real_list.pop_back();
    ideal_true_list[0].swap(ideal_true_list[1]);
    
    if(OnSameLine(real_list)){
        vector<vector<double> >().swap(RotMatix);
        vector<double>().swap(TranVect);

        vector<vector<double> >().swap(ideal_true_list);
        vector<vector<double> >().swap(real_list);
        return 0;
    }

    RotateCoor(ideal_true_list,real_list,RotMatix,TranVect);//real list atoms in a line
    if(atom_name==" C2'"&&residue.atoms[find_atom_in_residue(residue," O2'")].movable) ChangeCoor(tmp,RotMatix,TranVect,residue.atoms[find_atom_in_residue(residue," O2'")].xyz);
    if(atom_name==" C4'"&&residue.atoms[find_atom_in_residue(residue," C4'")].movable) ChangeCoor(tmp,RotMatix,TranVect,residue.atoms[find_atom_in_residue(residue," C4'")].xyz);

    if(atom_name==" C3'"&&residue.atoms[find_atom_in_residue(residue," C3'")].movable) ChangeCoor(tmp,RotMatix,TranVect,residue.atoms[find_atom_in_residue(residue," C3'")].xyz);
    if(atom_name==" C1'"&&residue.atoms[find_atom_in_residue(residue," C1'")].movable) ChangeCoor(tmp,RotMatix,TranVect,residue.atoms[find_atom_in_residue(residue," C1'")].xyz);
    vector<vector<double> >().swap(RotMatix);
    vector<double>().swap(TranVect);

    vector<vector<double> >().swap(ideal_true_list);
    vector<vector<double> >().swap(real_list);
    return 1;
}//Algorithm variant?

int check_P_chirality(ResidueUnit &residue_n,ResidueUnit &residue_n1){
    //cout<<"check_P_chirality"<<endl;
    vector<double> tmp(3,0);
    vector<vector<double> > ideal_true_list;
    vector<vector<double> > real_list;
    ideal_true_list.assign(4,tmp);
    real_list.assign(4,tmp);
    for(int i=0;i<4;i++){
        for(int j=0;j<3;j++){
            ideal_true_list[i][j]=P_chirality[i][j];
        }
    }

    for(int i=0;i<3;i++) real_list[0][i]=residue_n1.atoms[find_atom_in_residue(residue_n1," OP1")].xyz[i];
    for(int i=0;i<3;i++) real_list[1][i]=residue_n1.atoms[find_atom_in_residue(residue_n1," OP2")].xyz[i];
    for(int i=0;i<3;i++) real_list[2][i]=residue_n1.atoms[find_atom_in_residue(residue_n1," P  ")].xyz[i];
    for(int i=0;i<3;i++) real_list[3][i]=residue_n.atoms[find_atom_in_residue(residue_n," O3'")].xyz[i];

    vector<vector<double> > RotMatix;
    vector<double> TranVect;

    RotateCoor(ideal_true_list,real_list,RotMatix,TranVect);
    double SD_true=0;
    for(int i=0;i<4;i++){
        ChangeCoor(ideal_true_list[i],RotMatix,TranVect,tmp);
        SD_true+=Points2Distance2(real_list[i],tmp);
    }
    vector<vector<double> >().swap(RotMatix);
    vector<double>().swap(TranVect);

    ideal_true_list[0].swap(ideal_true_list[1]);
    RotateCoor(ideal_true_list,real_list,RotMatix,TranVect);
    double SD_false=0;
    for(int i=0;i<4;i++){
        ChangeCoor(ideal_true_list[i],RotMatix,TranVect,tmp);
        SD_false+=Points2Distance2(real_list[i],tmp);
    }
    vector<vector<double> >().swap(RotMatix);
    vector<double>().swap(TranVect);

    if(SD_true<SD_false){
        vector<vector<double> >().swap(ideal_true_list);
        vector<vector<double> >().swap(real_list);
        return 0;
    }
    
    vector<vector<double> >().swap(ideal_true_list);
    vector<vector<double> >().swap(real_list);
    if(residue_n1.atoms[find_atom_in_residue(residue_n1," OP1")].movable==0||residue_n1.atoms[find_atom_in_residue(residue_n1," OP1")].movable==0) return 0;//?
    residue_n1.atoms[find_atom_in_residue(residue_n1," OP1")].xyz.swap(residue_n1.atoms[find_atom_in_residue(residue_n1," OP2")].xyz);
    return 1;
}

int fix_chirality(ModelUnit &pep){
    int moved = 0;
    for (int c=0; c<pep.chains.size(); c++){
        for (int r=0; r<pep.chains[c].residues.size(); r++){
            if(pep.chains[c].residues[r].atoms.size()>=12){
                moved+=check_carbon_chirality(pep.chains[c].residues[r]," C2'");
                moved+=check_carbon_chirality(pep.chains[c].residues[r]," C4'");
                //if(check_carbon_chirality(pep.chains[c].residues[r]," C3'")>0) cout<<pep.chains[c].residues[r].resi<<"\t"<<pep.chains[c].residues[r].resn<<endl;
                moved+=check_carbon_chirality(pep.chains[c].residues[r]," C3'");
                
            }

            if(pep.chains[c].residues[r].atoms.size()>=13){
                moved+=check_carbon_chirality(pep.chains[c].residues[r]," C1'");
            }

            if(r != (pep.chains[c].residues.size()-1) && isConnected(pep,c,r) && pep.chains[c].residues[r].atoms.size()>=8&&pep.chains[c].residues[r+1].atoms.size()>=8){
                moved+=check_P_chirality(pep.chains[c].residues[r],pep.chains[c].residues[r+1]);
                if(check_P_chirality(pep.chains[c].residues[r],pep.chains[c].residues[r+1])>0) cout<<pep.chains[c].residues[r].resi<<"\t"<<pep.chains[c].residues[r].resn<<endl;
            }
        }
    }
    //if(moved > 0)cout<<moved<<" SwapChirality"<<endl;
    return moved;
    //return 0;
}

#endif