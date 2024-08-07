/* Purpose: Fix atomic clashes */

#ifndef AtomicClashes_HPP
#define AtomicClashes_HPP 1

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

// The order of each array is the standard PDB order, with hydrogen atoms omitted. VDW radii are in units of Angstroms.

// A
// Atomic order: P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N9 C8 N7 C5 C6 N6 N1 C2 N3 C4
double A_vdw[22] = {1.8, 1.52, 1.52, 1.52, 1.7, 1.7, 1.52, 1.7, 1.52, 1.7, 1.52, 1.7, 1.55, 1.7, 1.55, 1.7, 1.7, 1.55, 1.55, 1.7, 1.55, 1.7};

// C
// Atomic order: P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1 C2 O2 N3 C4 N4 C5 C6
double C_vdw[20] = {1.8, 1.52, 1.52, 1.52, 1.7, 1.7, 1.52, 1.7, 1.52, 1.7, 1.52, 1.7, 1.55, 1.7, 1.52, 1.55, 1.7, 1.55, 1.7, 1.7};

// G
// Atomic order: P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N9 C8 N7 C5 C6 O6 N1 C2 N2 N3 C4
double G_vdw[23] = {1.8, 1.52, 1.52, 1.52, 1.7, 1.7, 1.52, 1.7, 1.52, 1.7, 1.52, 1.7, 1.55, 1.7, 1.55, 1.7, 1.7, 1.52, 1.55, 1.7, 1.55, 1.55, 1.7};

// U
// Atomic order: P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1 C2 O2 N3 C4 O4 C5 C6
double U_vdw[20] = {1.8, 1.52, 1.52, 1.52, 1.7, 1.7, 1.52, 1.7, 1.52, 1.7, 1.52, 1.7, 1.55, 1.7, 1.52, 1.55, 1.7, 1.52, 1.7, 1.7};


void check_moveable(ModelUnit &pep, vector<vector<int> >&moveable_mat)
{

    vector<vector<int> >().swap(moveable_mat);
    size_t c,r,a;
    for (size_t c=0; c<pep.chains.size(); c++)
    {
        vector<int> tmp_vec(pep.chains[c].residues.size(),0);
        for (size_t r=0; r<pep.chains[c].residues.size(); r++)
        {
            for (size_t a=0; a<pep.chains[c].residues[r].atoms.size(); a++)
            {
                tmp_vec[r]+=pep.chains[c].residues[r].atoms[a].movable>0;
            }
        }
        moveable_mat.push_back(tmp_vec);
        vector<int>().swap(tmp_vec);
    }
    return;
}

int fix_clashes (ModelUnit &pep, vector<vector<int> >&moveable_mat){
	
	// Initialize variables
	double vdw1 = 0;
	double vdw2 = 0;

	int moved = 0;

	// OSP
	// Atomic order: P OP1 OP2 O5' C5' C4' O4' C3' O3'
	double OSP_bonds[9][4] = {
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{2.648, 0, 0, 0},
		{1.610, 2.505, 2.505, 2.504}
	};
     map<string,int> N;
    
    N[" P  "]=0; N[" OP1"]=1; N[" OP2"]=2; N[" O5'"]=3; N[" C5'"]=4; N[" C4'"]=5;
    N[" O4'"]=6; N[" C3'"]=7; N[" O3'"]=8; N[" C2'"]=9; N[" O2'"]=10; N[" C1'"]=11;
 	// Iterate through chains
	for (int c1=0; c1<pep.chains.size(); c1++){
		for (int c2=c1; c2<pep.chains.size(); c2++){	
			// Iterate through residues
			for (int r1=0; r1<pep.chains[c1].residues.size(); r1++){
				string resname1 = pep.chains[c1].residues[r1].resn;
				int resnumber1 = pep.chains[c1].residues[r1].resi;
				for (int r2=0; r2<pep.chains[c2].residues.size(); r2++){
                    if (moveable_mat[c1][r1] + moveable_mat[c2][r2]==0) continue;
					string resname2 = pep.chains[c2].residues[r2].resn;
					int resnumber2 = pep.chains[c2].residues[r2].resi;

					// Check each combination of r1 and r2 only once
					if (c1 == c2 && r1 >= r2) continue;


					/* Check C1'-C1' distance */
					if(pep.chains[c1].residues[r1].atoms.size()>=12&&pep.chains[c2].residues[r2].atoms.size()>=12){
					// Get (x,y,z) coordinates
					double C1x1 = pep.chains[c1].residues[r1].atoms[11].xyz[0];
					double C1y1 = pep.chains[c1].residues[r1].atoms[11].xyz[1];
					double C1z1 = pep.chains[c1].residues[r1].atoms[11].xyz[2];
					double C1x2 = pep.chains[c2].residues[r2].atoms[11].xyz[0];
					double C1y2 = pep.chains[c2].residues[r2].atoms[11].xyz[1];
					double C1z2 = pep.chains[c2].residues[r2].atoms[11].xyz[2];
					// Calculate the distance between two C1' atoms
					double C1_distance = sqrt((C1x2-C1x1)*(C1x2-C1x1) + (C1y2-C1y1)*(C1y2-C1y1) + (C1z2-C1z1)*(C1z2-C1z1));
					//cout << C1_distance << endl;

					if (C1_distance > 18.161) continue; // Not possible for a clash to occur
					}


						// Iterate through all pairs of atoms
						for (int a1=0; a1<pep.chains[c1].residues[r1].atoms.size(); a1++){
							string atom1 = pep.chains[c1].residues[r1].atoms[a1].name;
							int i1;
							if(pep.chains[c1].residues[r1].atoms.size()<12) i1=N[atom1];
							else i1=a1;
							for (int a2=0; a2<pep.chains[c2].residues[r2].atoms.size(); a2++){
								string atom2 = pep.chains[c2].residues[r2].atoms[a2].name;
								int i2;
								if(pep.chains[c2].residues[r2].atoms.size()<12) i2=N[atom2];
								else i2=a2;
								// Don't check OSP bonds
								if (r1+1==r2 && c1==c2 && i1 < 9 && i2 < 4 && OSP_bonds[i1][i2] > 0) continue;
									// Get the VDW radius of each atom in the pair
				 					if (resname1 == "  A"){
				 						vdw1 = A_vdw[i1];
				 					} else if (resname1 == "  C"){
				 						vdw1 = C_vdw[i1];
				 					} else if (resname1 == "  G"){
				 						vdw1 = G_vdw[i1];
				 					} else if (resname1 == "  U"){
				 						vdw1 = U_vdw[i1];
				 					}
				 					if (resname2 == "  A"){
				 						vdw2 = A_vdw[i2];
				 					} else if (resname2 == "  C"){
				 						vdw2 = C_vdw[i2];
				 					} else if (resname2 == "  G"){
				 						vdw2 = G_vdw[i2];
				 					} else if (resname2 == "  U"){
				 						vdw2 = U_vdw[i2];
				 					}
				 					// Calculate the clash distance
				 					double c = vdw1 + vdw2 - 0.4;
				 					// Get (x,y,z) coordinates
									double a1x = pep.chains[c1].residues[r1].atoms[a1].xyz[0];
									double a1y = pep.chains[c1].residues[r1].atoms[a1].xyz[1];
									double a1z = pep.chains[c1].residues[r1].atoms[a1].xyz[2];
									double a2x = pep.chains[c2].residues[r2].atoms[a2].xyz[0];
									double a2y = pep.chains[c2].residues[r2].atoms[a2].xyz[1];
									double a2z = pep.chains[c2].residues[r2].atoms[a2].xyz[2];
									// Calculate the atomic distance
									double d = sqrt((a2x-a1x)*(a2x-a1x) + (a2y-a1y)*(a2y-a1y) + (a2z-a1z)*(a2z-a1z));
					
									// Compare the atomic and clash distances
				 					if (d < c && d > 0.0001){
				 						//cout << "clash" << resname1 << resnumber1 << atom1 << resname2 << resnumber2 << atom2 << "distance= " << d << ' ' << "vdw= " << c << endl;
				 						// Calculate difference from VDW distance
				 						double delta = d - (vdw1 + vdw2); 
				 						if (pep.chains[c1].residues[r1].atoms[a1].movable==1 && pep.chains[c2].residues[r2].atoms[a2].movable==0){
											moved++;
											pep.chains[c1].residues[r1].atoms[a1].xyz[0] = (delta/d)*(a2x-a1x) + a1x;
											pep.chains[c1].residues[r1].atoms[a1].xyz[1] = (delta/d)*(a2y-a1y) + a1y;
											pep.chains[c1].residues[r1].atoms[a1].xyz[2] = (delta/d)*(a2z-a1z) + a1z;
										} else if (pep.chains[c1].residues[r1].atoms[a1].movable==0 && pep.chains[c2].residues[r2].atoms[a2].movable==1){
											moved++;
											pep.chains[c2].residues[r2].atoms[a2].xyz[0] = (delta/d)*(a1x-a2x) + a2x;
											pep.chains[c2].residues[r2].atoms[a2].xyz[1] = (delta/d)*(a1y-a2y) + a2y;
											pep.chains[c2].residues[r2].atoms[a2].xyz[2] = (delta/d)*(a1z-a2z) + a2z;
										} else if (pep.chains[c1].residues[r1].atoms[a1].movable==1 && pep.chains[c2].residues[r2].atoms[a2].movable==1){
											moved++;
											// move each atom half the required distance
											pep.chains[c2].residues[r2].atoms[a2].xyz[0] = 0.5*((delta/d)*(a1x-a2x)) + a2x;
											pep.chains[c2].residues[r2].atoms[a2].xyz[1] = 0.5*((delta/d)*(a1y-a2y)) + a2y;
											pep.chains[c2].residues[r2].atoms[a2].xyz[2] = 0.5*((delta/d)*(a1z-a2z)) + a2z;
											pep.chains[c1].residues[r1].atoms[a1].xyz[0] = 0.5*((delta/d)*(a2x-a1x)) + a1x;
											pep.chains[c1].residues[r1].atoms[a1].xyz[1] = 0.5*((delta/d)*(a2y-a1y)) + a1y;
											pep.chains[c1].residues[r1].atoms[a1].xyz[2] = 0.5*((delta/d)*(a2z-a1z)) + a1z;
										}
										//cout << "clash" << resname1 << resnumber1 << atom1 << resname2 << resnumber2 << atom2 << "distance= " << d << ' ' << "vdw= " << c << endl;
				 					}
				 			}
				 		}
 				}
 			}
 		}
 	}
 	return moved;
}

//6.1

vector<ResidueUnit> find_clashes(ModelUnit &pep,vector<vector<ResidueUnit> > &new_resis){
	// Initialize variables
	double vdw1 = 0;
	double vdw2 = 0;

	size_t clash=0;
	size_t least=0;
	int best=0;
	// OSP
	// Atomic order: P OP1 OP2 O5' C5' C4' O4' C3' O3'
	double OSP_bonds[9][4] = {
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{2.648, 0, 0, 0},
		{1.610, 2.505, 2.505, 2.504}
	};
     map<string,int> N;
    
    N[" P  "]=0; N[" OP1"]=1; N[" OP2"]=2; N[" O5'"]=3; N[" C5'"]=4; N[" C4'"]=5;
    N[" O4'"]=6; N[" C3'"]=7; N[" O3'"]=8; N[" C2'"]=9; N[" O2'"]=10; N[" C1'"]=11;
    
	cout<<new_resis[0][0].resi<<endl;
	for(int i=0;i<new_resis.size();i++){
		for(size_t j=0;j<new_resis[i].size();j++){
			for(int c=0;c<pep.chains.size();c++){
				for(size_t r=0;r<pep.chains[c].residues.size();r++){
					string resname1 = new_resis[i][j].resn;
					int resnumber1 = new_resis[i][j].resi;
					string resname2 = pep.chains[c].residues[r].resn;
				    int resnumber2 = pep.chains[c].residues[r].resi;
					for (int a1=0; a1<new_resis[i][j].atoms.size(); a1++){
						string atom1 = new_resis[i][j].atoms[a1].name;
						int i1;
						cout<<atom1<<endl;
						if(new_resis[i][j].atoms.size()<12) i1=N[atom1];
						else i1=a1;
						for(int a2=0; a2<pep.chains[c].residues[r].atoms.size(); a2++){
							string atom2 = pep.chains[c].residues[r].atoms[a2].name;
							int i2=a2;

							        if (resname1 == "  A"){
				 						vdw1 = A_vdw[i1];
				 					} else if (resname1 == "  C"){
				 						vdw1 = C_vdw[i1];
				 					} else if (resname1 == "  G"){
				 						vdw1 = G_vdw[i1];
				 					} else if (resname1 == "  U"){
				 						vdw1 = U_vdw[i1];
				 					}
				 					if (resname2 == "  A"){
				 						vdw2 = A_vdw[i2];
				 					} else if (resname2 == "  C"){
				 						vdw2 = C_vdw[i2];
				 					} else if (resname2 == "  G"){
				 						vdw2 = G_vdw[i2];
				 					} else if (resname2 == "  U"){
				 						vdw2 = U_vdw[i2];
				 					}
									// Calculate the clash distance
				 					double ca = vdw1 + vdw2 - 0.4;
				 					// Get (x,y,z) coordinates
									double a1x =new_resis[i][j].atoms[a1].xyz[0];
									double a1y =new_resis[i][j].atoms[a1].xyz[1];
									double a1z =new_resis[i][j].atoms[a1].xyz[2];
									double a2x = pep.chains[c].residues[r].atoms[a2].xyz[0];
									double a2y = pep.chains[c].residues[r].atoms[a2].xyz[1];
									double a2z = pep.chains[c].residues[r].atoms[a2].xyz[2];
									// Calculate the atomic distance
									double d = sqrt((a2x-a1x)*(a2x-a1x) + (a2y-a1y)*(a2y-a1y) + (a2z-a1z)*(a2z-a1z));
					
									// Compare the atomic and clash distances
				 					if (d < ca && d > 0.0001) clash++;
									if(atom1==" C4'"&&atom2=="C4'"){
										if(d<6.5) clash++;
								
							}
						}
					}
				}
			}
		}
		if(i==0) least=clash;
		else if(clash<=least){
			least=clash;
			best=i;
		}
		clash=0;
	}
	vector<ResidueUnit> new_resi;
	for(size_t o=0;o<new_resis[best].size();o++){
		new_resi.push_back(new_resis[best][o]);
	}
	return new_resi;
}

vector<int> atom_clash_C4(ModelUnit &pep,vector<vector<ResidueUnit> > &new_resis){
	int clash=0;
	vector<int> clashes;
	for(int h=0;h<new_resis.size();h++){
		clash=0;
		for(int i=0;i<new_resis[h].size();i++){
			for(int j=0;j<pep.chains.size();j++){
				for(size_t k=0;k<pep.chains[j].residues.size();k++){
					if(Points2Distance(new_resis[h][i].atoms[0].xyz,pep.chains[j].residues[k].atoms[5].xyz)<6.1) clash++;
				}
			}
		}
		clashes.push_back(clash);
	}
	return clashes;
}

vector<ResidueUnit> atom_clashes_C4_quick(ModelUnit &pep, vector<vector<ResidueUnit> > &new_resis, ResidueUnit resi_f, ResidueUnit resi_r){
    double r=6;
    vector<AtomUnit> atoms_concerned;
    for(int i=0;i<pep.chains.size();i++){
        for(size_t j=0;j<pep.chains[i].residues.size();j++){
            for(int k=0;k<pep.chains[i].residues[j].atoms.size();k++){
                if(Points2Distance(resi_f.atoms[5].xyz,pep.chains[i].residues[j].atoms[k].xyz)
                +Points2Distance(resi_r.atoms[5].xyz,pep.chains[i].residues[j].atoms[k].xyz)
                <6.2*(resi_r.resi-resi_f.resi)+2*r)
                atoms_concerned.push_back(pep.chains[i].residues[j].atoms[k]);
            }
        }
    }
    /*int clash;
    vector<int> clashes;
    for(size_t i=0;i<atoms_concerned.size();i++){
        for(int j=0;j<new_resis.size();j++){
            clash=0;
            for(int k=0;k<new_resis[j].size();k++){
                if(Points2Distance(atoms_concerned[i].xyz,new_resis[j][k].atoms[0].xyz)<r) clash++; 
            }
            clashes.push_back(clash);
        }
    }
    vector<AtomUnit>().swap(atoms_concerned);
    return clashes;*/
	int clash;
    size_t least=atoms_concerned.size();
	int best=0;
    for(size_t i=0;i<atoms_concerned.size();i++){
        for(int j=0;j<new_resis.size();j++){
            clash=0;
            for(int k=0;k<new_resis[j].size();k++){
                if(Points2Distance(atoms_concerned[i].xyz,new_resis[j][k].atoms[0].xyz)<r) clash++; 
            }
			if(least>clash){
				least=clash;
				best=j;
			}
        }
    }
    vector<AtomUnit>().swap(atoms_concerned);
	return new_resis[best];
}
#endif
