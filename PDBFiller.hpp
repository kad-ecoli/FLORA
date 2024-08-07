/* Purpose1: find gaps in PDB structures and fill them through multiple methods. */
#ifndef PDBFiller_HPP
#define PDBFiller_HPP 1

#include <cmath>
#include <random>
#include <ctime>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include "PDBParser.hpp"
#include "GeometryTools.hpp"
#include "Superpose.hpp"
#include "AtomicClashes.hpp"
#include "MissingRNAatom.hpp"
#include "AdjustPosition.hpp"

using namespace std;

typedef vector<pair<size_t,bool> > sizet_bool_vc;
typedef vector<pair<size_t,size_t> > sizet_sizet_vc;



struct Model_F{
    vector<string> chains_F;
};

//dev random(user random call)
static unsigned int X=123456789, Y=987654321, Z=43219876, C=6543217;

unsigned int devrand(void){
    int fn;
    unsigned int r;

    fn=open("/dev/urandom",O_RDONLY);
    if(fn==-1) exit(-1);
    if(read(fn,&r,4)!=4) exit(-1);
    close(fn);

    return r;
}

void init_KISS(){
    X=devrand();
    while(!(Y=devrand()));
    Z=devrand();

    C=devrand()%698769068+1;
}

unsigned int JKISS()
{
    unsigned long long t;

    X=314527869*X+1234567;
    Y^=Y<<5; Y^=Y>>7;Y^=Y<<22;
    t=429458439ULL*Z+C; C=t>>32;Z=t;

    return X+Y+Z;
}

Model_F read_fasta_structure(const char* filename){
    Model_F model;
    string filename_str=(string)filename;
    ifstream fp;
    fp.open(filename,ios::in);
    string line;
    int chainI=-1;
    while(fp.good()){
        getline(fp,line);
        if(!fp.good()) break;
        if(line[0]=='>'){
            chainI++;
            continue;
        }
        if(model.chains_F.size()==chainI) model.chains_F.push_back(line);
        else if(model.chains_F.size()>chainI) model.chains_F[chainI]+=line;
    }
    fp.close();
    if(chainI==-1) cout<<"No valid sequence in the file!"<<endl;
    return model;
}

void print_fasta_structure(Model_F &model){
	for(size_t i=0;i<model.chains_F.size();i++){
		cout<<">:"<<i<<endl;
		cout<<model.chains_F[i]<<endl;
	}
}

int compare_fasta_pdb(const string sequence, ChainUnit &chain_p, int chainI){
	int unpaired_res=0;
	int paired_res=0;
	int unknown_res=0;
	size_t i,j;
	for(i=0;i<sequence.size();i++){
		for(j=0;j<chain_p.residues.size();j++){
			if(chain_p.residues[j].het==1) continue;
			if(i==chain_p.residues[j].resi-1){
				if(tolower(sequence[i])!=aa3to1(chain_p.residues[j].resn,2)){
					unpaired_res++;
				}else paired_res++;
			}
		}
	}
	unknown_res=chain_p.residues.size()-unpaired_res-paired_res;
	if(unpaired_res) cout<<"Error: the PDB file cannot match the sequence. ("<<unpaired_res<<" residues in chain "<<chainI<<")\n";
	if(unknown_res) cout<<"Error: "<<unknown_res<<" residues are not included in chain "<<chainI<<".\n";
	return unpaired_res+unknown_res;
}

int compare_fasta_pdb(Model_F &model_f, ModelUnit &model_p){
	if(model_f.chains_F.size()<model_p.chains.size()){
        cerr<<"Error: some chains are not included in the sequence file!\n";
        return 1;
    }
    else if(model_f.chains_F.size()>model_p.chains.size()){
        //cout<<model_p.chains.size()<<endl;
        cerr<<"Error: missing chains in the PDB file\n";
        return 1;
    }
    int e=0;
    for(int i=0;i<model_f.chains_F.size();i++){
        e+=compare_fasta_pdb(model_f.chains_F[i],model_p.chains[i],i+1);
	}
	return e;
}

void write_SEQRES(const char* outfile, Model_F &fasta, ModelUnit &pdb){
    if(!compare_fasta_pdb(fasta,pdb)){
        ofstream fp(outfile);
        for(int i=0;i<fasta.chains_F.size();i++){
            int num_line=fasta.chains_F[i].size()/13+1;
            for(int j=0;j<num_line;j++){
                fp<<"SEQRES"<<" "
                  <<setw(3)<<j+1<<" "
                  <<pdb.chains[i].chainID<<" "
                  <<setw(4)<<fasta.chains_F[i].size()<<" ";
                for(int k=0;13*j+k<fasta.chains_F[i].size()&&k<13;k++){
                    fp<<setw(4)<<toupper(fasta.chains_F[i][13*j+k]);
                }
                fp<<endl;
            }
        }
        fp.close();
    }
}

void modify_RNA_structure(ModelUnit &pdb_in){
    for(int c=0;c<pdb_in.chains.size();c++){
        for(int r=0;r<pdb_in.chains[c].residues.size();r++){
            if(pdb_in.chains[c].residues[r].icode=='1') pdb_in.chains[c].residues[r].icode=' ';
        }
    }
}

inline int find_atom(ResidueUnit &resi,string atom_name){
    for(int i=0;i<resi.atoms.size();i++){
        if(resi.atoms[i].name==atom_name) return i;
    }
    return -1;
}

sizet_sizet_vc gapinfo(const string sequence,ChainUnit &chain_p, const int option=0){
    //option=0, include the gap at both ends; option=1, only record inner gaps
    sizet_sizet_vc gaps;
    size_t start;
    if(option==0) start=1;
    if(option==1) start=chain_p.residues[0].resi;
    for(size_t i=0;i<chain_p.residues.size();i++){
        if(chain_p.residues[i].resi==start){
            start++;
        }else{
            gaps.push_back(pair<size_t,size_t>(start,chain_p.residues[i].resi));
            start=chain_p.residues[i].resi+1;
        }
    }
    if(option==0&&start<sequence.size()+1) gaps.push_back(pair<size_t,size_t>(start,sequence.size()+1));
    return gaps;
}

sizet_sizet_vc gapinfo(const string sequence,ChainUnit &chain_p, vector<vector<ResidueUnit> > &gapif,const int option=0,const int thre=4){
    sizet_sizet_vc gaps=gapinfo(sequence,chain_p,option);
    if(option==1||option==0){
        //gap.first is the first one in the gap, but gap.second is (the last one in the gap + 1)
        for(int g=0;g<gaps.size();g++){
            size_t f,r;
            if(g==0) f=gaps[g].first-1;
            else f=gaps[g].first-gaps[g-1].second;
            if(g==gaps.size()-1) r=sequence.size()-gaps[g].second+1;
            else r=gaps[g+1].first-gaps[g].second;
            cout<<f<<"\t"<<r<<endl;
            if(f>thre) f=thre;
            if(r>thre) r=thre;
            bool a=false;
            bool b=false;
            vector<ResidueUnit> res_f;//e.g. 5 is missing, get 1,2,3,4
            vector<ResidueUnit> res_r;//e.g. 9 is missing, get 13,12,11,10
            for(size_t i=0;i<chain_p.residues.size();i++){
                if(!a&&f>0&&chain_p.residues[i].resi==gaps[g].first-f){
                    for(int j=0;j<f;j++){
                        res_f.push_back(chain_p.residues[i+j]);
                    }
                    a=true;
                }
                if(!b&&r>0&&chain_p.residues[i].resi==gaps[g].second-1+r){
                    for(int j=0;j<r;j++){
                        res_r.push_back(chain_p.residues[i-j]);
                    }
                    b=true;
                }
                if(a&&b) break;
            }
            gapif.push_back(res_f);
            gapif.push_back(res_r);
            vector<ResidueUnit>().swap(res_f);
            vector<ResidueUnit>().swap(res_r);
        }
    }
    return gaps;
}

int fill_residue(const string sequence, ChainUnit &chain_p){
    size_t i,j;
    sizet_sizet_vc gaps=gapinfo(sequence,chain_p,0);
    /*for(i=0;i<gaps.size();i++){
        cout<<chain_p.chainID_full<<gaps[i].first<<" "<<gaps[i].second<<endl;
    }*/
    string resn_full="   ";
    ResidueUnit Residue;
    Residue.het=false;
    Residue.icode='0';
    Residue.resn="   ";
    for(i=0;i<gaps.size();i++){
        for(j=gaps[i].first;j<gaps[i].second;j++){
            Residue.resi=j;
            Residue.resn[2]=toupper(sequence[j-1]);
            chain_p.residues.push_back(Residue);
        }
    }
    return 1;
}

void fill_residue(Model_F &model_f,ModelUnit &model_p){
    standardize_residue_order(model_p);
    for(size_t i=0;i<model_f.chains_F.size();i++){
        fill_residue(model_f.chains_F[i],model_p.chains[i]);
    }
    standardize_residue_order(model_p);
}

void fill_new_residue(ChainUnit &chain, vector<ResidueUnit> &new_resi){
    for(size_t i=0;i<new_resi.size();i++){
        chain.residues.push_back(new_resi[i]);
    }
    standardize_residue_order(chain);
}

void replace_new_residue(ChainUnit &chain, vector<ResidueUnit> &new_resi){
    int pointer=0;
    standardize_residue_order(new_resi);
    for(size_t i=0;i<new_resi.size();i++){
        for(size_t j=pointer;j<chain.residues.size();j++){
            if(chain.residues[j].resi==new_resi[i].resi){
                pointer=j;
                chain.residues[j].het=new_resi[i].het;
                chain.residues[j].icode=new_resi[i].icode;
                chain.residues[j].resn=new_resi[i].resn;
                vector<AtomUnit>().swap(chain.residues[j].atoms);
                for(int k=0;k<new_resi[i].atoms.size();k++) chain.residues[j].atoms.push_back(new_resi[i].atoms[k]);
                break;
            }
        }
    }
    standardize_residue_order(chain);
}

void delete_moveable_residue(ChainUnit &chain){
    int moveable_num=0;
    for(int r=0;r<chain.residues.size();r++){
        int confidence=0;
        for(int a=0;a<chain.residues[r].atoms.size();a++){
            if(chain.residues[r].atoms[a].movable<1) confidence++;
        }
        if(!confidence&&chain.residues[r].icode!='1'){
            chain.residues[r].icode='0';
            moveable_num++;
        }
    }
    
    for(int r=0;r<chain.residues.size();r++){
        if (chain.residues[r].icode!='0') continue;
        
        for (int j=r+1;j<chain.residues.size();j++)
        {
            if (chain.residues[j].icode=='0') continue;
            chain.residues[r].het  =chain.residues[j].het;
            chain.residues[r].resi =chain.residues[j].resi;
            chain.residues[r].icode=chain.residues[j].icode;
            chain.residues[r].resn =chain.residues[j].resn;
            chain.residues[r].atoms=chain.residues[j].atoms;
            chain.residues[j].icode='0';
            break;
        }
    }
    for(int num=0;num<moveable_num;num++) chain.residues.pop_back();
}

bool full_PDB(ModelUnit &pdb_in,Model_F &fasta_in){
    if(pdb_in.chains.size()!=fasta_in.chains_F.size()) return false;
    for(int c=0;c<pdb_in.chains.size();c++){
        for(int r=1;r<pdb_in.chains[c].residues.size();r++){
            if(pdb_in.chains[c].residues[r].resi-pdb_in.chains[c].residues[r-1].resi!=1) return false;
        }
        string fasta_o=pdb2fasta(pdb_in.chains[c]);
        for(int j=0;j<fasta_o.size();j++) fasta_o[j]=tolower(fasta_o[j]);
        string fasta_i=fasta_in.chains_F[c];
        for(int j=0;j<fasta_i.size();j++) fasta_i[j]=tolower(fasta_i[j]);
        if(strcmp(fasta_o.c_str(),fasta_i.c_str())) return false;
    }
    return true;
}

bool fill_linear(string sequence, ResidueUnit resi_f, ResidueUnit resi_r, vector<ResidueUnit> &new_resi){
    ResidueUnit resi_n;
    resi_n.het=false;
    resi_n.icode=' ';
    resi_n.resn="   ";
    AtomUnit atom;
    atom.movable=1;
    atom.name="    ";
    atom.xyz.assign(3,0);
    resi_n.atoms.assign(12,atom);
    if(resi_r.resi<resi_f.resi){
        cout<<"invalid input!"<<endl;
        return false;
    }
    for(int i=1;i<resi_r.resi-resi_f.resi;i++){
        resi_n.resi=resi_f.resi+i;
        resi_n.resn[2]=toupper(sequence[resi_n.resi-1]);
        for(int j=0;j<12;j++){
            resi_n.atoms[j].name=resi_f.atoms[j].name;
            //cout<<resi_n.atoms[j].name;
            for(int k=0;k<3;k++){
                resi_n.atoms[j].xyz[k]=
                resi_f.atoms[j].xyz[k]+
                (1.0*i*(resi_r.atoms[j].xyz[k]-resi_f.atoms[j].xyz[k]))/(resi_r.resi-resi_f.resi);
                //cout<<" "<<resi_n.atoms[j].xyz[k];
            }
            //cout<<endl;
        }  
        new_resi.push_back(resi_n);
    }
    
    return true;
}

int fill_linear(const string sequence, ChainUnit &chain_p){
    vector<vector<ResidueUnit> > gapif;
    vector<ResidueUnit> new_resi;
    
    gapinfo(sequence,chain_p,gapif,1,1);

    for(size_t i=0;i<gapif.size();i+=2){
        if(gapif[i].size()>0&&gapif[i+1].size()>0){
            fill_linear(sequence,gapif[i][0],gapif[i+1][0],new_resi);
            fill_new_residue(chain_p,new_resi);
            vector<ResidueUnit>().swap(new_resi);
        }
    }
    vector<vector<ResidueUnit> >().swap(gapif);
    return 1;
}


inline double y_l(double d, int N, double h){
    int h_2_d4=4*h*h+d*d;
    double theta=(M_PI-acos(2*h/sqrt(h_2_d4)))/N;
    double y=sqrt(h_2_d4)*sin(theta);
    return y;
}

inline double y1_l(double d, int N, double h){
    int x_2_d4=4*h*h+d*d;
    int theta=(acos(2*h/sqrt(x_2_d4))-M_PI)/N;
    double y1=sqrt(x_2_d4)*(2*d*cos(theta)-4*N*h*sin(theta))/(N*x_2_d4);
    return y1;
}

inline double x_ykb(double y, double k, double b){
    return (y-b)/k;
}

inline double x_ykxy(double y, double k, double x0, double y0){
    return ((y-y0)/k+x0);
}

void RNAArc(double distance, int gap_size, int start_resi, string sequence, vector<vector<double> > &adjacent_atoms, vector<ResidueUnit> &new_resi){
    //gap_size=atom_between+1
    double intercept=distance*sin(M_PI_2/gap_size);
    double k_a=2*sin(M_PI/gap_size);
    double h;
    if(intercept>6.2){
        h=x_ykxy(5.8,y1_l(distance,gap_size,0),0,y_l(distance,gap_size,0));
        //cout<<distance/gap_size<<" "<<h<<" "<<y1_l(distance,gap_size,0)<<" "<<y_l(distance,gap_size,0)<<endl;
        while(y_l(distance,gap_size,h)>6.2){
            h=x_ykxy(5.8,y1_l(distance,gap_size,h),h,y_l(distance,gap_size,h));
            //cout<<"1 "<<h<<endl;
        }
    }
    else if(intercept<5.8){
        h=x_ykb(5.8,k_a,0);
        //cout<<distance/gap_size<<" "<<h<<" "<<y1_l(distance,gap_size,0)<<" "<<y_l(distance,gap_size,0)<<endl;
        while(y_l(distance,gap_size,h)>6.2){
            h=x_ykxy(5.8,y1_l(distance,gap_size,h),h,y_l(distance,gap_size,h));
            //cout<<"2 "<<h<<endl;
        }
    }
    else h=0;
    //cout<<"A"<<endl;
    ResidueUnit resi_temp;
    AtomUnit atom_temp0;
    AtomUnit atom_temp1;
    AtomUnit atom_temp2;
    vector<double> endpoint_m;
    vector<double> endpoint_p;
    vector<double> center;
    vector<double> loc;
    vector<double> coor;
    vector<double> coor_1;
    vector<double> coor_2;
    endpoint_m.assign(3,0);endpoint_m[0]-=distance/2;
    endpoint_p.assign(3,0);endpoint_p[0]+=distance/2;
    center.assign(3,0);center[1]=h;
    loc.assign(3,0);loc[1]=-1;
    coor.assign(3,0);
    coor_1.assign(3,0);
    coor_2.assign(3,0);
    vector<vector<double> > coors;
    double R=Points2Distance(center,endpoint_p);
    double phi=acos(h/R)-M_PI_2;
    double theta=(2*M_PI-2*acos(h/R))/gap_size;
    for(int i=1;i<gap_size;i++){
        double angle=phi+i*theta;
        coor[0]=R*cos(angle);
        coor[1]=R*sin(angle)+h;
        coors.push_back(coor);
        coor_1[0]=coor[0]+1.03*cos(angle)-1.10*sin(angle);
        coor_1[1]=coor[1]+1.03*sin(angle)+1.10*cos(angle);
        coors.push_back(coor_1);
        coor_2[0]=coor[0]+1.03*cos(angle)+2.14*sin(angle);
        coor_2[1]=coor[1]+1.03*sin(angle)-2.14*cos(angle);
        coors.push_back(coor_2);
    }

    vector<double> mid1;
    vector<double> mid2;
    mid1.assign(3,0);
    mid2.assign(3,0);
    for(int i=0;i<3;i++) mid1[i]=(adjacent_atoms[1][i]+adjacent_atoms[2][i])/2;
    for(int i=0;i<3;i++) mid2[i]=(adjacent_atoms[0][i]+adjacent_atoms[3][i])/2;
    vector<double> gap;
    gap.assign(3,0);
    vector<double> h_vector;
    h_vector.assign(3,0);
    subtract(adjacent_atoms[1],adjacent_atoms[2],gap);
    subtract(mid1,mid2,h_vector);
    //ShowMyvector(gap);
    //ShowMyvector(h_vector);
    //cout<<innerproduct(gap,h_vector)<<"\t"<<innerproduct(gap,gap)<<endl;
    multi(innerproduct(gap,h_vector)/innerproduct(gap,gap),gap);
    //ShowMyvector(gap);
    subtract(h_vector,gap,h_vector);
    //ShowMyvector(h_vector);
    norm_warnless(h_vector);
    //ShowMyvector(h_vector);
    for(int i=0;i<3;i++) mid1[i]-=h_vector[i];

    vector<vector<double> > xyz_list1;
    vector<vector<double> > xyz_list2;
    vector<vector<double> > RotMatix;
    vector<double> TranVect;
    xyz_list1.push_back(endpoint_p);
    xyz_list1.push_back(endpoint_m);
    xyz_list1.push_back(loc);
    xyz_list2.push_back(adjacent_atoms[1]);
    xyz_list2.push_back(adjacent_atoms[2]);
    xyz_list2.push_back(mid1);
    RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);
    resi_temp.het=false;
    resi_temp.icode=' ';
    resi_temp.resn="   ";
    atom_temp0.movable=1;
    atom_temp0.name=" C4'";
    atom_temp0.xyz.assign(3,0);
    atom_temp1.movable=1;
    atom_temp1.name=" C5'";
    atom_temp1.xyz.assign(3,0);
    atom_temp2.movable=1;
    atom_temp2.name=" O3'";
    atom_temp2.xyz.assign(3,0);
    for(int i=0;i<coors.size();i+=3){
        resi_temp.resi=start_resi+1+i/3;
        resi_temp.resn[2]=toupper(sequence[resi_temp.resi-1]);
        ChangeCoor(coors[i],RotMatix,TranVect,atom_temp0.xyz);
        ChangeCoor(coors[i+1],RotMatix,TranVect,atom_temp1.xyz);
        ChangeCoor(coors[i+2],RotMatix,TranVect,atom_temp2.xyz);
        resi_temp.atoms.push_back(atom_temp1);
        resi_temp.atoms.push_back(atom_temp0);
        resi_temp.atoms.push_back(atom_temp2);
        new_resi.push_back(resi_temp);
        vector<AtomUnit>().swap(resi_temp.atoms);
    }

    vector<double>().swap(endpoint_m);
    vector<double>().swap(endpoint_p);
    vector<double>().swap(center);
    vector<double>().swap(loc);
    vector<double>().swap(coor);
    vector<double>().swap(coor_1);
    vector<double>().swap(coor_2);
    vector<vector<double> >().swap(coors);
    vector<double>().swap(mid1);
    vector<double>().swap(mid2);
    vector<double>().swap(h_vector);
    vector<double>().swap(gap);
    vector<vector<double> >().swap(xyz_list1);
    vector<vector<double> >().swap(xyz_list2);
    vector<vector<double> >().swap(RotMatix);
    vector<double>().swap(TranVect);
}

bool fill_arc(const string sequence, ChainUnit &chain_p){
    vector<vector<ResidueUnit> > gapif;
    vector<ResidueUnit> new_resi;
    
    gapinfo(sequence,chain_p,gapif,0,2);
    cout<<gapif.size()/2<<endl;
    cout<<gapif[0].size()<<"\t"<<gapif[1].size()<<endl;
    for(size_t i=0;i<gapif.size();i+=2){
        //cout<<i<<endl;
        if(gapif[i].size()==0||gapif[i+1].size()==0) continue;
        if(gapif[i].size()==2&&gapif[i+1].size()==2&&Points2Distance(gapif[i][1].atoms[find_atom(gapif[i][1]," C4'")].xyz,gapif[i+1][1].atoms[find_atom(gapif[i+1][1]," C4'")].xyz)<5.4*(gapif[i+1][1].resi-gapif[i][1].resi)){
            vector<vector<double> > adjacent_atoms;
            adjacent_atoms.push_back(gapif[i][0].atoms[find_atom(gapif[i][0]," C4'")].xyz);
            adjacent_atoms.push_back(gapif[i][1].atoms[find_atom(gapif[i][1]," C4'")].xyz);
            adjacent_atoms.push_back(gapif[i+1][1].atoms[find_atom(gapif[i+1][1]," C4'")].xyz);
            adjacent_atoms.push_back(gapif[i+1][0].atoms[find_atom(gapif[i+1][0]," C4'")].xyz);
            //cout<<gapif[i][1].resi<<endl;
            RNAArc(Points2Distance(adjacent_atoms[1],adjacent_atoms[2]),gapif[i+1][1].resi-gapif[i][1].resi,gapif[i][1].resi,sequence,adjacent_atoms,new_resi);
            vector<vector<double> >().swap(adjacent_atoms);
        }
        else{
            fill_linear(sequence,gapif[i][gapif[i].size()-1],gapif[i+1][gapif[i+1].size()-1],new_resi);
        }
        fill_new_residue(chain_p,new_resi);
        vector<ResidueUnit>().swap(new_resi);
    }
    vector<vector<ResidueUnit> >().swap(gapif);
    return 1;
}

/*FLORA_6 - arc filling*/
void fill_arc(Model_F &model_f,ModelUnit &model_p){
    bool i1;
    for(size_t i=0;i<model_p.chains.size();i++){
        //i1=fill_SRS(model_f.chains_F[i],model_p.chains[i]);
        i1=fill_arc(model_f.chains_F[i],model_p.chains[i]);
    }
    standardize_residue_order(model_p);
}


void extrapolate_both1(const string sequence, ChainUnit &chain_p,vector<ResidueUnit> &ref_list_p,vector<ResidueUnit> &ref_list_n,vector<ResidueUnit> &new_residue){
    int p=ref_list_p.size();
    int n=ref_list_n.size();
    pair<size_t,size_t> gap;
    gap.first=ref_list_p[ref_list_p.size()-1].resi+1;
    gap.second=ref_list_n[ref_list_n.size()-1].resi;
    vector<double> tmp(3,0);
    vector<vector<double> > xyz_listr;
    vector<vector<double> > xyz_listt;
    vector<vector<double> > RotMatix;
    vector<double> TranVect;
    ResidueUnit Residue;
    Residue.het=false;
    Residue.icode=' ';
    AtomUnit Atom;
    Atom.movable=1;
    Atom.xyz.assign(3,0);
    Residue.atoms.assign(12,Atom);
    string resn_full="   ";
    xyz_listr.assign(36,tmp);
    xyz_listt.assign(36,tmp);
    for(size_t i=0;i<gap.second-gap.first;i++){
        if(i%2){
            for(int j=0;j<3;j++){
                for(int k=0;k<12;k++){
                    xyz_listt[12*j+k][0]=ref_list_n[n-2-j].atoms[k].xyz[0];
                    xyz_listt[12*j+k][1]=ref_list_n[n-2-j].atoms[k].xyz[1];
                    xyz_listt[12*j+k][2]=ref_list_n[n-2-j].atoms[k].xyz[2];
                    xyz_listr[12*j+k][0]=ref_list_n[n-1-j].atoms[k].xyz[0];
                    xyz_listr[12*j+k][1]=ref_list_n[n-1-j].atoms[k].xyz[1];
                    xyz_listr[12*j+k][2]=ref_list_n[n-1-j].atoms[k].xyz[2];
                }
            }
            RotateCoor(xyz_listt,xyz_listr, RotMatix, TranVect);
            Residue.resi=ref_list_n[n-1].resi-1;
            resn_full[2]=toupper(sequence[Residue.resi-1]);
            Residue.resn=resn_full;
            for(int k=0;k<12;k++){
                Residue.atoms[k].name=ref_list_n[n-1].atoms[k].name;
                ChangeCoor(ref_list_n[n-1].atoms[k].xyz,RotMatix,TranVect,Residue.atoms[k].xyz);
            }
            ref_list_n.push_back(Residue);
            n++;
            new_residue.push_back(Residue);
        }else{
            for(int j=0;j<3;j++){
                for(int k=0;k<12;k++){
                    xyz_listt[12*j+k][0]=ref_list_p[p-2-j].atoms[k].xyz[0];
                    xyz_listt[12*j+k][1]=ref_list_p[p-2-j].atoms[k].xyz[1];
                    xyz_listt[12*j+k][2]=ref_list_p[p-2-j].atoms[k].xyz[2];
                    xyz_listr[12*j+k][0]=ref_list_p[p-1-j].atoms[k].xyz[0];
                    xyz_listr[12*j+k][1]=ref_list_p[p-1-j].atoms[k].xyz[1];
                    xyz_listr[12*j+k][2]=ref_list_p[p-1-j].atoms[k].xyz[2];
                }
            }
            RotateCoor(xyz_listt,xyz_listr, RotMatix, TranVect);
            Residue.resi=ref_list_p[p-1].resi+1;
            resn_full[2]=toupper(sequence[Residue.resi-1]);
            Residue.resn=resn_full;
            for(int k=0;k<12;k++){
                Residue.atoms[k].name=ref_list_p[p-1].atoms[k].name;
                ChangeCoor(ref_list_p[p-1].atoms[k].xyz,RotMatix,TranVect,Residue.atoms[k].xyz);
            }
            ref_list_p.push_back(Residue);
            p++;
            new_residue.push_back(Residue);
        }
    }
    vector<ResidueUnit>().swap(ref_list_p);
    vector<ResidueUnit>().swap(ref_list_n);
    vector<vector<double> >().swap(xyz_listt);
    vector<vector<double> >().swap(xyz_listr);
}

void extrapolate_prev(const string sequence, ChainUnit &chain_p,vector<ResidueUnit> &ref_list_p,int length,vector<ResidueUnit> &new_residue){
    pair<size_t,size_t> gap;
    gap.first=ref_list_p[ref_list_p.size()-1].resi+1;
    gap.second=gap.first+length;
    int p=ref_list_p.size();
    vector<double> tmp(3,0);
    vector<vector<double> > xyz_listr;
    vector<vector<double> > xyz_listt;
    vector<vector<double> > RotMatix;
    vector<double> TranVect;
    ResidueUnit Residue;
    Residue.het=false;
    Residue.icode=' ';
    AtomUnit Atom;
    Atom.movable=1;
    Atom.xyz.assign(3,0);
    Residue.atoms.assign(12,Atom);
    string resn_full="   ";
    xyz_listr.assign(36,tmp);
    xyz_listt.assign(36,tmp);
    for(size_t i=0;i<gap.second-gap.first;i++){
            for(int j=0;j<3;j++){
                for(int k=0;k<12;k++){
                    xyz_listt[12*j+k][0]=ref_list_p[p-2-j].atoms[k].xyz[0];
                    xyz_listt[12*j+k][1]=ref_list_p[p-2-j].atoms[k].xyz[1];
                    xyz_listt[12*j+k][2]=ref_list_p[p-2-j].atoms[k].xyz[2];
                    xyz_listr[12*j+k][0]=ref_list_p[p-1-j].atoms[k].xyz[0];
                    xyz_listr[12*j+k][1]=ref_list_p[p-1-j].atoms[k].xyz[1];
                    xyz_listr[12*j+k][2]=ref_list_p[p-1-j].atoms[k].xyz[2];
                }
            }
            RotateCoor(xyz_listt,xyz_listr, RotMatix, TranVect);
            Residue.resi=ref_list_p[p-1].resi+1;
            resn_full[2]=toupper(sequence[Residue.resi-1]);
            Residue.resn=resn_full;
            for(int k=0;k<12;k++){
                Residue.atoms[k].name=ref_list_p[p-1].atoms[k].name;
                ChangeCoor(ref_list_p[p-1].atoms[k].xyz,RotMatix,TranVect,Residue.atoms[k].xyz);
            }
            ref_list_p.push_back(Residue);
            p++;
            new_residue.push_back(Residue);
    }
}

void extrapolate_next(const string sequence, ChainUnit &chain_p,vector<ResidueUnit> &ref_list_n,int length,vector<ResidueUnit> &new_residue){
    int n=ref_list_n.size();
    pair<size_t,size_t> gap;
    gap.second=ref_list_n[ref_list_n.size()-1].resi;
    gap.first=gap.second-length;
    vector<double> tmp(3,0);
    vector<vector<double> > xyz_listr;
    vector<vector<double> > xyz_listt;
    vector<vector<double> > RotMatix;
    vector<double> TranVect;
    ResidueUnit Residue;
    Residue.het=false;
    Residue.icode=' ';
    AtomUnit Atom;
    Atom.movable=1;
    Atom.xyz.assign(3,0);
    Residue.atoms.assign(12,Atom);
    string resn_full="   ";
    xyz_listr.assign(36,tmp);
    xyz_listt.assign(36,tmp);
    for(size_t i=0;i<gap.second-gap.first;i++){
            for(int j=0;j<3;j++){
                for(int k=0;k<12;k++){
                    xyz_listt[12*j+k][0]=ref_list_n[n-2-j].atoms[k].xyz[0];
                    xyz_listt[12*j+k][1]=ref_list_n[n-2-j].atoms[k].xyz[1];
                    xyz_listt[12*j+k][2]=ref_list_n[n-2-j].atoms[k].xyz[2];
                    xyz_listr[12*j+k][0]=ref_list_n[n-1-j].atoms[k].xyz[0];
                    xyz_listr[12*j+k][1]=ref_list_n[n-1-j].atoms[k].xyz[1];
                    xyz_listr[12*j+k][2]=ref_list_n[n-1-j].atoms[k].xyz[2];
                }
            }
            RotateCoor(xyz_listt,xyz_listr, RotMatix, TranVect);
            Residue.resi=ref_list_n[n-1].resi-1;
            resn_full[2]=toupper(sequence[Residue.resi-1]);
            Residue.resn=resn_full;
            for(int k=0;k<12;k++){
                Residue.atoms[k].name=ref_list_n[n-1].atoms[k].name;
                ChangeCoor(ref_list_n[n-1].atoms[k].xyz,RotMatix,TranVect,Residue.atoms[k].xyz);
            }
            ref_list_n.push_back(Residue);
            n++;
            new_residue.push_back(Residue);
    }
}

int extrapolate(const string sequence, ChainUnit &chain_p){
    vector<vector<ResidueUnit> > gapif;
    sizet_sizet_vc gaps=gapinfo(sequence,chain_p,gapif,0,4);
    cout<<chain_p.residues.size()<<endl;
    vector<ResidueUnit> new_resi;
    
    int i1=0;
    int length;
    for(size_t g=0;g<gapif.size();g+=2){
        length=gaps[g/2].second-gaps[g/2].first;
        if(gapif[g].size()>=4&&gapif[g+1].size()>=4){
            //cout<<"extrapolate from both sides"<<endl;
            extrapolate_both1(sequence,chain_p,gapif[g],gapif[g+1],new_resi);
            //cout<<chain_p.residues.size()<<endl;
        }
        else if(gapif[g].size()>=4&&gapif[g+1].size()<4){
            //cout<<"extrapolate from the front"<<endl;
            extrapolate_prev(sequence,chain_p,gapif[g],length,new_resi);
            //cout<<length<<endl;
            //cout<<chain_p.residues.size()<<endl;
        }
        else if(gapif[g].size()<4&&gapif[g+1].size()>=4){
            //cout<<"extrapolate from the rear"<<endl;
            extrapolate_next(sequence,chain_p,gapif[g+1],length,new_resi);
            //cout<<length<<endl;
            //cout<<chain_p.residues.size()<<endl;
        }
        else i1++;
        //cout<<new_resi.size()<<endl;
        //cout<<"gap No."<<g<<endl;
    }
    //if(i1==0){
    for(size_t i=0;i<new_resi.size();i++){
        chain_p.residues.push_back(new_resi[i]);
    }
    //}
    
    vector<ResidueUnit>().swap(new_resi);
    standardize_residue_order(chain_p);
    cout<<chain_p.residues.size()<<endl;
    return i1;
}

/*FLORA_1 - extrapolation*///Every residue in the model should have atoms on the backbone
void extrapolate(Model_F &model_f,ModelUnit &model_p){
    
    for(size_t i=0;i<model_f.chains_F.size();i++){
        int i1=3;
        while(i1){
            extrapolate(model_f.chains_F[i],model_p.chains[i]);
            i1--;
        }
        
    }
    standardize_residue_order(model_p);
}

inline bool too_far(double distance,int atom_bet){
    double mean;
    double SD;
    /*mean=9.8687*log(1.0+atom_bet)+5.0694;
    SD=0.7508*(1.0+atom_bet)-0.217;*///0
    /*mean=9.8188*log(1.0+atom_bet)+4.83;
    SD=0.742*(1.0+atom_bet)-0.242;*///3
    /*mean=9.9969*log(1.0+atom_bet)+5.4004;
    SD=0.769*(1.0+atom_bet)-0.2284;*///4
    /*mean=9.9111*log(1.0+atom_bet)+5.3799;
    SD=0.8048*(1.0+atom_bet)-0.3608;*///5
    /*mean=9.7385*log(1.0+atom_bet)+4.9579;
    SD=0.7925*(1.0+atom_bet)-0.4516;*///7
    /*mean=9.9634*log(1.0+atom_bet)+5.422;
    SD=0.7752*(1.0+atom_bet)-0.2626;*///8
    if(distance/(atom_bet+1)>4.5&&distance/(atom_bet+1)<5.2) return false;
    //if(distance<mean+3*SD&&distance>mean-3*SD) return false;
    else return true;
}

inline double dis_fest(int atom_bet){
    double mean;
    double SD;
    mean=9.7385*log(0.0+atom_bet);
    SD=0.7925*(0.0+atom_bet)-0.4516;
    return (6.5*atom_bet);
}

bool fill_gap_new(const string sequence,pair<size_t,size_t> gap,ChainUnit &chain_p,int prev,int next,vector<ResidueUnit> &new_resi,const int option=5){
    vector<ResidueUnit> ref_list_p;
    vector<ResidueUnit> ref_list_n;
    if(prev>4) prev=4;
    if(next>4) next=4;
    bool a=false;
    bool b=false;
    for(size_t i=0;i<chain_p.residues.size();i++){
        if(!a && chain_p.residues[i].resi==gap.first-prev){
            for(int j=0;j<prev;j++){
                ref_list_p.push_back(chain_p.residues[i+j]);
            }
            a=true;
        }
        if(!b && chain_p.residues[i].resi==gap.second-1+next){
            for(int j=0;j<next;j++){
                ref_list_n.push_back(chain_p.residues[i-j]);
            }
            b=true;
        }
        if(a&&b) break;
    }
    int p=ref_list_p.size();
    int n=ref_list_n.size();
    if(p+n<2) return false;
    vector<double> tmp(3,0);
    vector<vector<double> > xyz_listr;
    vector<vector<double> > xyz_listt;
    vector<vector<double> > RotMatix;
    vector<double> TranVect;
    ResidueUnit Residue_p;
    ResidueUnit Residue_n;
    ResidueUnit Residue_l;
    vector<ResidueUnit> Residues_p;
    vector<ResidueUnit> Residues_n;
    Residue_p.het=false;
    Residue_n.het=false;
    Residue_l.het=false;
    Residue_p.icode=' ';
    Residue_n.icode=' ';
    Residue_l.icode=' ';
    AtomUnit Atom;
    Atom.movable=0;
    Atom.name="  P ";
    Atom.xyz.assign(3,0);
    Residue_p.atoms.assign(12,Atom);
    Residue_n.atoms.assign(12,Atom);
    Residue_l.atoms.assign(12,Atom);
    string resn_full="   ";
    if(p>0) Residues_p.push_back(ref_list_p[ref_list_p.size()-1]);
    if(n>0) Residues_n.push_back(ref_list_n[ref_list_n.size()-1]);
    if(p>1){
        xyz_listr.assign(12*(p-1),tmp);
        xyz_listt.assign(12*(p-1),tmp);
        for(int j=0;j<p-1;j++){
            for(int k=0;k<12;k++){
                xyz_listt[12*j+k][0]=ref_list_p[p-2-j].atoms[k].xyz[0];
                xyz_listt[12*j+k][1]=ref_list_p[p-2-j].atoms[k].xyz[1];
                xyz_listt[12*j+k][2]=ref_list_p[p-2-j].atoms[k].xyz[2];
                xyz_listr[12*j+k][0]=ref_list_p[p-1-j].atoms[k].xyz[0];
                xyz_listr[12*j+k][1]=ref_list_p[p-1-j].atoms[k].xyz[1];
                xyz_listr[12*j+k][2]=ref_list_p[p-1-j].atoms[k].xyz[2];
            }
        }
        RotateCoor(xyz_listt,xyz_listr, RotMatix, TranVect);
        for(size_t i=gap.first;i<gap.second;i++){
            Residue_p.resi=i;
            resn_full[2]=toupper(sequence[Residue_p.resi-1]);
            Residue_p.resn=resn_full;
            for(int j=0;j<12;j++){
                Residue_p.atoms[j].movable=1;
                Residue_p.atoms[j].name=ref_list_p[0].atoms[j].name;
                ChangeCoor(Residues_p[Residues_p.size()-1].atoms[j].xyz,RotMatix,TranVect,Residue_p.atoms[j].xyz);
            }
            Residues_p.push_back(Residue_p);
        }
        vector<vector<double> >().swap(xyz_listt);
        vector<vector<double> >().swap(xyz_listr);
    }
    if(n>1){
        xyz_listr.assign(12*(n-1),tmp);
        xyz_listt.assign(12*(n-1),tmp);
        for(int j=0;j<n-1;j++){
            for(int k=0;k<12;k++){
                xyz_listt[12*j+k][0]=ref_list_n[n-2-j].atoms[k].xyz[0];
                xyz_listt[12*j+k][1]=ref_list_n[n-2-j].atoms[k].xyz[1];
                xyz_listt[12*j+k][2]=ref_list_n[n-2-j].atoms[k].xyz[2];
                xyz_listr[12*j+k][0]=ref_list_n[n-1-j].atoms[k].xyz[0];
                xyz_listr[12*j+k][1]=ref_list_n[n-1-j].atoms[k].xyz[1];
                xyz_listr[12*j+k][2]=ref_list_n[n-1-j].atoms[k].xyz[2];
            }
        }
        RotateCoor(xyz_listt,xyz_listr, RotMatix, TranVect);
        for(size_t i=gap.second-1;i>=gap.first;i--){
            Residue_n.resi=i;
            resn_full[2]=toupper(sequence[Residue_n.resi-1]);
            Residue_n.resn=resn_full;
            for(int j=0;j<12;j++){
                Residue_n.atoms[j].movable=1;
                Residue_n.atoms[j].name=ref_list_n[0].atoms[j].name;
                ChangeCoor(Residues_n[Residues_n.size()-1].atoms[j].xyz,RotMatix,TranVect,Residue_n.atoms[j].xyz);
            }
            Residues_n.push_back(Residue_n);
        }
        vector<vector<double> >().swap(xyz_listt);
        vector<vector<double> >().swap(xyz_listr);
    }
    if(p==0) for(size_t j=1;j<Residues_n.size();j++) new_resi.push_back(Residues_n[j]);
    if(n==0) for(size_t j=1;j<Residues_p.size();j++) new_resi.push_back(Residues_p[j]);
    if(p==0||n==0){
        vector<ResidueUnit>().swap(Residues_p);
        vector<ResidueUnit>().swap(Residues_n);
        vector<ResidueUnit>().swap(ref_list_p);
        vector<ResidueUnit>().swap(ref_list_n);
        return true;
    } 
    int atom_bet=0;
    double distance;
    vector<pair<size_t,size_t> > ij_pair;
    for(size_t i=0;i<Residues_p.size();i++){
        for(size_t j=0;j<Residues_n.size();j++){
            atom_bet=gap.second-gap.first-i-j;
            if(option==4&&atom_bet!=0) continue;//when option=4, only chains that meet could be recorded.
            else if(atom_bet<0) continue;
            distance=Points2Distance(Residues_p[i].atoms[0].xyz,Residues_n[j].atoms[0].xyz);
            if(!too_far(distance,atom_bet)) ij_pair.push_back(pair<size_t,size_t>(i,j));
        }
    }
    cout<<ij_pair.size()<<endl;
    pair<size_t,size_t> gap_l;
    if(ij_pair.size()>0&&option!=2){
        size_t pointer=0;
        atom_bet=gap.second-gap.first;
        for(size_t i=0;i<ij_pair.size();i++){
            if(atom_bet>=gap.second-gap.first-ij_pair[i].first-ij_pair[i].second){
                atom_bet=gap.second-gap.first-ij_pair[i].first-ij_pair[i].second;
                pointer=i;
            }
        }
        gap_l.first=ij_pair[pointer].first;
        gap_l.second=ij_pair[pointer].second;
    }else{
        gap_l.first=0;
        gap_l.second=0;
    }
    for(size_t j=1;j<=gap_l.first;j++) new_resi.push_back(Residues_p[j]);
    for(size_t j=1;j<=gap_l.second;j++) new_resi.push_back(Residues_n[j]);
    double begin,end;
    for(size_t i=1;i<=gap.second-gap.first-gap_l.second-gap_l.first;i++){
        Residue_l.resi=Residues_p[gap_l.first].resi+i;
        resn_full[2]=toupper(sequence[Residue_l.resi-1]);
        Residue_l.resn=resn_full;
        for(int j=0;j<12;j++){
            Residue_l.atoms[j].movable=1;
            Residue_l.atoms[j].name=Residues_p[gap_l.first].atoms[j].name;
            for(int k=0;k<3;k++){
                begin=Residues_p[gap_l.first].atoms[j].xyz[k];
                end=Residues_n[gap_l.second].atoms[j].xyz[k];
                Residue_l.atoms[j].xyz[k]=begin+i*(end-begin)/(gap.second-gap.first-gap_l.second-gap_l.first+1);
            }
        }
        new_resi.push_back(Residue_l);
    }
    vector<pair<size_t,size_t> >().swap(ij_pair);
    vector<ResidueUnit>().swap(Residues_p);
    vector<ResidueUnit>().swap(Residues_n);
    vector<ResidueUnit>().swap(ref_list_p);
    vector<ResidueUnit>().swap(ref_list_n);
    vector<vector<double> >().swap(xyz_listt);
    vector<vector<double> >().swap(xyz_listr);
    return true;
}

int fill_gap_new(const string sequence,ChainUnit &chain_p,const int option=5){
    sizet_sizet_vc gaps;
    size_t start=1;
    for(size_t i=0;i<chain_p.residues.size();i++){
        if(chain_p.residues[i].resi==start){
            start++;
        }else{
            gaps.push_back(pair<size_t,size_t>(start,chain_p.residues[i].resi));
            start=chain_p.residues[i].resi+1;
        }
    }
    if(start<sequence.size()+1) gaps.push_back(pair<size_t,size_t>(start,sequence.size()+1));
    vector<ResidueUnit> new_resi;
    sizet_sizet_vc gaps_temp;
    int gaps_size=gaps.size();
    while(gaps_size){
        for(size_t g=0;g<gaps.size();g++){
            int ep,en;
            bool finish;
            if(g==0) ep=gaps[g].first-1;
            else ep=gaps[g].first-gaps[g-1].second;
            if(g==gaps.size()-1) en=sequence.size()+1-gaps[g].second;
            else en=gaps[g+1].first-gaps[g].second;
            finish=fill_gap_new(sequence,gaps[g],chain_p,ep,en,new_resi,option);
            if(!finish) gaps_temp.push_back(gaps[g]);
            else{
                for(size_t j=0;j<new_resi.size();j++){
                    chain_p.residues.push_back(new_resi[j]);
                }
                standardize_residue_order(chain_p);
                vector<ResidueUnit>().swap(new_resi);
            }
        }
        gaps.swap(gaps_temp);
        sizet_sizet_vc().swap(gaps_temp);
        if(gaps_size==gaps.size()) break;
        gaps_size=gaps.size();
    }
    return 1;
}

/*FLORA_2 - Combination Method I  */
/*FLORA_4 - Combination Method III*/
/*FLORA_5 - Combination Method IV */
void fill_gap_new(Model_F &model_f,ModelUnit &model_p,const int option=5){
    int i1;
    for(size_t i=0;i<model_p.chains.size();i++){
        i1=fill_gap_new(model_f.chains_F[i],model_p.chains[i],option);
    }
    standardize_residue_order(model_p);
}

bool fill_gap(const string sequence,pair<size_t,size_t> gap,ChainUnit &chain_p,int prev,int next,vector<ResidueUnit> &new_resi){
    vector<ResidueUnit> ref_list_p;
    vector<ResidueUnit> ref_list_n;
    if(prev>4) prev=4;
    if(next>4) next=4;
    bool a=false;
    bool b=false;
    for(size_t i=0;i<chain_p.residues.size();i++){
        if(!a && chain_p.residues[i].resi==gap.first-prev){
            for(int j=0;j<prev;j++){
                ref_list_p.push_back(chain_p.residues[i+j]);
            }
            a=true;
        }
        if(!b && chain_p.residues[i].resi==gap.second-1+next){
            for(int j=0;j<next;j++){
                ref_list_n.push_back(chain_p.residues[i-j]);
            }
            b=true;
        }
        if(a&&b) break;
    }
    int p=ref_list_p.size();
    int n=ref_list_n.size();
    if(p+n<2) return false;
    vector<double> tmp(3,0);
    vector<vector<double> > xyz_listr;
    vector<vector<double> > xyz_listt;
    vector<vector<double> > RotMatix;
    vector<double> TranVect;
    ResidueUnit Residue_p;
    ResidueUnit Residue_n;
    Residue_p.icode='0';
    Residue_n.icode='0';
    AtomUnit Atom;
    Atom.xyz.assign(3,0);
    string resn_full="   ";
    int option=0;
    int thre_p=0;
    int thre_n=0;
    for(size_t i=0;i<gap.second-gap.first;i++){
        double dis_p,dis_n=0;
        if(p>1){
            (p>4)?thre_p=4:thre_p=p;
            xyz_listr.assign(12*(thre_p-1),tmp);
            xyz_listt.assign(12*(thre_p-1),tmp);
            for(int j=0;j<thre_p-1;j++){
                for(int k=0;k<12;k++){
                    xyz_listt[12*j+k][0]=ref_list_p[p-2-j].atoms[k].xyz[0];
                    xyz_listt[12*j+k][1]=ref_list_p[p-2-j].atoms[k].xyz[1];
                    xyz_listt[12*j+k][2]=ref_list_p[p-2-j].atoms[k].xyz[2];
                    xyz_listr[12*j+k][0]=ref_list_p[p-1-j].atoms[k].xyz[0];
                    xyz_listr[12*j+k][1]=ref_list_p[p-1-j].atoms[k].xyz[1];
                    xyz_listr[12*j+k][2]=ref_list_p[p-1-j].atoms[k].xyz[2];
                }
            }
            RotateCoor(xyz_listt,xyz_listr, RotMatix, TranVect);
            Residue_p.het=false;
            Residue_p.resi=ref_list_p[p-1].resi+1;
            resn_full[2]=toupper(sequence[Residue_p.resi-1]);
            Residue_p.resn=resn_full;
            Residue_p.icode=' ';
            Atom.movable=1;
            for(int k=0;k<12;k++){
                Atom.name=ref_list_p[p-1].atoms[k].name;
                ChangeCoor(ref_list_p[p-1].atoms[k].xyz,RotMatix,TranVect,Atom.xyz);
                Residue_p.atoms.push_back(Atom);
            }
            if(n>0){
                dis_p=Points2Distance(Residue_p.atoms[7].xyz,ref_list_n[n-1].atoms[7].xyz);
                if(dis_p>dis_fest(ref_list_n[n-1].resi-Residue_p.resi)){
                    Residue_p.icode='0';
                    Residue_p.atoms.clear();
                    Residue_p.atoms.shrink_to_fit();
                }
                if(dis_p>Points2Distance(ref_list_p[p-1].atoms[7].xyz,ref_list_n[n-1].atoms[7].xyz)){
                    Residue_p.icode='0';
                    Residue_p.atoms.clear();
                    Residue_p.atoms.shrink_to_fit();
                }
            }
        }
        if(n>1){
            (n>4)?thre_n=4:thre_n=n;
            xyz_listr.assign(12*(thre_n-1),tmp);
            xyz_listt.assign(12*(thre_n-1),tmp);
            for(int j=0;j<thre_n-1;j++){
                for(int k=0;k<12;k++){
                    xyz_listt[12*j+k][0]=ref_list_n[n-2-j].atoms[k].xyz[0];
                    xyz_listt[12*j+k][1]=ref_list_n[n-2-j].atoms[k].xyz[1];
                    xyz_listt[12*j+k][2]=ref_list_n[n-2-j].atoms[k].xyz[2];
                    xyz_listr[12*j+k][0]=ref_list_n[n-1-j].atoms[k].xyz[0];
                    xyz_listr[12*j+k][1]=ref_list_n[n-1-j].atoms[k].xyz[1];
                    xyz_listr[12*j+k][2]=ref_list_n[n-1-j].atoms[k].xyz[2];
                }
            }
            RotateCoor(xyz_listt,xyz_listr, RotMatix, TranVect);
            Residue_n.het=false;
            Residue_n.resi=ref_list_n[n-1].resi-1;
            resn_full[2]=toupper(sequence[Residue_n.resi-1]);
            Residue_n.resn=resn_full;
            Residue_n.icode=' ';
            Atom.movable=1;
            for(int k=0;k<12;k++){
                Atom.name=ref_list_n[n-1].atoms[k].name;
                ChangeCoor(ref_list_n[n-1].atoms[k].xyz,RotMatix,TranVect,Atom.xyz);
                Residue_n.atoms.push_back(Atom);
            }
            if(p>0){
                dis_n=Points2Distance(Residue_n.atoms[7].xyz,ref_list_p[p-1].atoms[7].xyz);
                if(dis_n>dis_fest((Residue_n.resi-ref_list_p[p-1].resi))){
                    Residue_n.icode='0';
                    Residue_n.atoms.clear();
                    Residue_n.atoms.shrink_to_fit();
                }
                if(dis_n>Points2Distance(ref_list_n[n-1].atoms[7].xyz,ref_list_p[p-1].atoms[7].xyz)){
                    Residue_n.icode='0';
                    Residue_n.atoms.clear();
                    Residue_n.atoms.shrink_to_fit();
                }
            }
        }
        if(Residue_p.icode!='0'&&Residue_n.icode!='0'){
            if(dis_p>dis_n){
                Residue_p.icode='0';
                Residue_p.atoms.clear();
                Residue_p.atoms.shrink_to_fit();
            }else{
                Residue_n.icode='0';
                Residue_n.atoms.clear();
                Residue_n.atoms.shrink_to_fit();
            }
        }
        if(Residue_p.icode==' '){
            ref_list_p.push_back(Residue_p);
            new_resi.push_back(Residue_p);
            Residue_p.atoms.clear();
            Residue_p.atoms.shrink_to_fit();
        }
        else if(Residue_n.icode==' '){
            ref_list_n.push_back(Residue_n);
            new_resi.push_back(Residue_n);
            Residue_n.atoms.clear();
            Residue_n.atoms.shrink_to_fit();
        }
        else{
            Residue_p.het=false;
            Residue_p.resi=ref_list_p[p-1].resi+1;
            resn_full[2]=toupper(sequence[Residue_p.resi-1]);
            Residue_p.resn=resn_full;
            Residue_p.icode=' ';
            Atom.movable=1;
            for(int k=0;k<12;k++){
                Atom.name=ref_list_p[p-1].atoms[k].name;
                double begin=0;
                double end=0;
                for(int l=0;l<3;l++){
                    begin=ref_list_p[p-1].atoms[k].xyz[l];
                    end=ref_list_n[n-1].atoms[k].xyz[l];
                    Atom.xyz[l]=begin+(end-begin)/(ref_list_n[n-1].resi-ref_list_p[p-1].resi);
                }
                Residue_p.atoms.push_back(Atom);
            }
            ref_list_p.push_back(Residue_p);
            new_resi.push_back(Residue_p);
        }
        vector<AtomUnit>().swap(Residue_p.atoms);
        vector<AtomUnit>().swap(Residue_n.atoms);
        p=ref_list_p.size();
        n=ref_list_n.size();
    }
    vector<ResidueUnit>().swap(ref_list_p);
    vector<ResidueUnit>().swap(ref_list_n);
    vector<vector<double> >().swap(xyz_listt);
    vector<vector<double> >().swap(xyz_listr);
    return true;
}

int fill_gap(const string sequence,ChainUnit &chain_p){
    sizet_sizet_vc gaps;
    size_t start=1;
    for(size_t i=0;i<chain_p.residues.size();i++){
        if(chain_p.residues[i].resi==start){
            start++;
        }else{
            gaps.push_back(pair<size_t,size_t>(start,chain_p.residues[i].resi));
            start=chain_p.residues[i].resi+1;
        }
    }
    if(start<sequence.size()+1) gaps.push_back(pair<size_t,size_t>(start,sequence.size()+1));
    vector<ResidueUnit> new_resi;
    sizet_sizet_vc gaps_temp;
    int gaps_size=gaps.size();
    while(gaps_size){
        for(size_t g=0;g<gaps.size();g++){
            int ep,en;
            bool finish;
            if(g==0) ep=gaps[g].first-1;
            else ep=gaps[g].first-gaps[g-1].second;
            if(g==gaps.size()-1) en=sequence.size()+1-gaps[g].second;
            else en=gaps[g+1].first-gaps[g].second;
            finish=fill_gap(sequence,gaps[g],chain_p,ep,en,new_resi);
            if(!finish) gaps_temp.push_back(gaps[g]);
            else{
                for(size_t j=0;j<new_resi.size();j++){
                    chain_p.residues.push_back(new_resi[j]);
                }
                standardize_residue_order(chain_p);
                vector<ResidueUnit>().swap(new_resi);
            }
        }
        gaps.swap(gaps_temp);
        sizet_sizet_vc().swap(gaps_temp);
        if(gaps_size==gaps.size()) break;
        gaps_size=gaps.size();
    }
    return 1;
}

/*FLORA_3 - Combination Method II*/
void fill_gap(Model_F &model_f,ModelUnit &model_p){
    int i1;
    for(size_t i=0;i<model_p.chains.size();i++){
        i1=fill_gap(model_f.chains_F[i],model_p.chains[i]);
    }
    standardize_residue_order(model_p);
}

int trim_backbone(ChainUnit &chain_p){
    /*string P=" P  ";string OP1=" OP1";string OP2=" OP2";string O5=" O5'";
    string C5=" C5'";string C4=" C4'";string O4=" O4'";string C3=" C3'";
    string O3=" O3'";string C2=" C2'";string O2=" O2'";string C1=" C1'";*/
    int omission=0;
    for(size_t i=0;i<chain_p.residues.size();i++){
        if(chain_p.residues[i].atoms.size()>12) continue;
        for(int j=0;j<12;j++){
            if(chain_p.residues[i].atoms[j-omission].name==" P  "||chain_p.residues[i].atoms[j-omission].name==" O5'"
             ||chain_p.residues[i].atoms[j-omission].name==" C5'"||chain_p.residues[i].atoms[j-omission].name==" C4'"
             ||chain_p.residues[i].atoms[j-omission].name==" C3'"||chain_p.residues[i].atoms[j-omission].name==" O3'") continue;
            for(int k=j-omission;k<11-omission;k++){
                chain_p.residues[i].atoms[k].name=chain_p.residues[i].atoms[k+1].name;
                chain_p.residues[i].atoms[k].xyz.swap(chain_p.residues[i].atoms[k+1].xyz);
            }
            omission++;
        }
        for(int l=0;l<omission;l++) chain_p.residues[i].atoms.pop_back();
        omission=0;
    }
    return 0;
}

void trim_backbone(ModelUnit &model_p){
    int i1=0;
    for(size_t i=0;i<model_p.chains.size();i++){
        i1=trim_backbone(model_p.chains[i]);
    }
    standardize_residue_order(model_p);
}

void SRS(vector<double> Sc,double r,vector<double> &xyz){
    uniform_real_distribution<double> interval(-1.0,1.0);
    normal_distribution<double> Gaussian(0.848,0.123);
    default_random_engine engine(JKISS());
    double costheta;
    double phi;
    do{costheta=Gaussian(engine);} while (costheta>1||costheta<-1);
    phi=M_PI+M_PI*interval(engine);
    xyz[0]=Sc[0]+r*sqrt(1-costheta*costheta)*cos(phi);
    xyz[1]=Sc[1]+r*sqrt(1-costheta*costheta)*sin(phi);
    xyz[2]=Sc[2]+r*costheta;
}

void new_coor(vector<double> &Sc,double r,double costheta,double phi,vector<double> &xyz){
    xyz[0]=Sc[0]+r*sqrt(1-costheta*costheta)*cos(phi);
    xyz[1]=Sc[1]+r*sqrt(1-costheta*costheta)*sin(phi);
    xyz[2]=Sc[2]+r*costheta;
}

void new_coor_Rot(vector<double> &Sc,double r,double costheta,double phi,vector<double> &xyz, vector<double> &RotVec){
    //cout<<RotVec[0]<<"\t"<<RotVec[1]<<"\t"<<RotVec[2]<<endl;
    double R2=RotVec[0]*RotVec[0]+RotVec[1]*RotVec[1]+RotVec[2]*RotVec[2];
    double XOZ=acos(RotVec[2]/sqrt(R2));
    double YOX=atan2(RotVec[1],RotVec[0]);
    cout<<phi<<"\t"<<costheta<<endl;
    double x=r*sqrt(1-costheta*costheta)*cos(phi);
    double y=r*sqrt(1-costheta*costheta)*sin(phi);
    double z=r*costheta;
    
    cout<<x<<"\t"<<y<<"\t"<<z<<endl;
    double xoz;
    if(z==0) xoz=0;
    else xoz=atan2(x,z);
    double rxoz=sqrt(x*x+z*z);
    x=rxoz*sin(xoz+XOZ);
    z=rxoz*cos(xoz+XOZ);
    
    double yox;
    if(x==0) yox=0;
    else yox=atan2(y,x);
    double ryox=sqrt(x*x+y*y);
    y=ryox*sin(yox+YOX);
    x=ryox*cos(yox+YOX);

    xyz[0]=Sc[0]+x;
    xyz[1]=Sc[1]+y;
    xyz[2]=Sc[2]+z;
}
//gapunit long stringent
int too_far_SRS(double distance,int gap_unit){
    double mean;
    double SD;
    mean=9.9111*log(gap_unit)+5.3799;
    SD=0.8048*(gap_unit)-0.3608;
    if(distance<mean+3*SD&&distance>mean-3*SD&&distance>5.2){
        //cout<<distance<<" "<<gap_unit<<endl;
        return 0;//
    } 
    else return 1;
}

bool ok_SRS(vector<ResidueUnit> &res_f,vector<ResidueUnit> &res_r,vector<double> atom_xyz,int i){
    int p=0;
    int resi;
    if(i%2){
        resi=res_f[0].resi+i/2+1;
    }else resi=res_r[0].resi-i/2;
    //cout<<i<<endl;
    /*cout<<resi<<endl;
    cout<<Points2Distance(res_f[0].atoms[0].xyz,res_r[0].atoms[0].xyz)<<endl;
    cout<<Points2Distance(res_f[0].atoms[0].xyz,atom_xyz)<<endl;
    cout<<Points2Distance(atom_xyz,res_r[0].atoms[0].xyz)<<endl;*/
    for(int j=0;j<res_f.size();j++){
        if(too_far_SRS(Points2Distance(res_f[j].atoms[0].xyz,atom_xyz),resi-res_f[j].resi)) return false;
    }
    for(int j=0;j<res_r.size();j++){
        if(too_far_SRS(Points2Distance(res_r[j].atoms[0].xyz,atom_xyz),res_r[j].resi-resi)) return false;
    }
    return true;
}

double inter(void){
    double x=0, y=0;
    x=JKISS()/4294967296.0;
    y=2*x-1;
    return y;
}

bool fill_SRS_terminal(string sequence, ChainUnit &chain_p,int f,int r,vector<ResidueUnit> &new_resi){
    vector<double> RotVec;
    RotVec.assign(3,0);
    RotVec[2]=1;
    ResidueUnit resi_new;
    resi_new.het=false;
    resi_new.icode=' ';
    resi_new.resn="   ";
    ResidueUnit resi_last;
    resi_last.het=false;
    resi_last.icode=' ';
    resi_last.resn="   ";
    AtomUnit atom;
    atom.movable=1;
    atom.name=" C4'";
    atom.xyz.assign(3,0);
    resi_new.atoms.assign(1,atom);
    resi_last.atoms.assign(1,atom);
    normal_distribution<double> Gaussian(0.848,0.123);
    uniform_real_distribution<double> interval(-1.0,1.0);
    default_random_engine engine;
    engine.seed(4922995424);
    double costheta;
    double phi;
    if(r==-1&&f!=-1){
        resi_new.resi=chain_p.residues[f].resi;
        resi_new.resn=chain_p.residues[f].resn;
        int f_C4=0;
        for(int j=0;j<chain_p.residues[f].atoms.size();j++){
            if(chain_p.residues[f].atoms[j].name==" C4'") f_C4=j;
        }
        resi_new.atoms[0].xyz[0]=chain_p.residues[f].atoms[f_C4].xyz[0];
        resi_new.atoms[0].xyz[1]=chain_p.residues[f].atoms[f_C4].xyz[1];
        resi_new.atoms[0].xyz[2]=chain_p.residues[f].atoms[f_C4].xyz[2];
        new_resi.push_back(resi_new);

        if(f!=0&&chain_p.residues[f].resi-chain_p.residues[f-1].resi==1){
            int f_1_C4=0;
            for(int l=0;l<chain_p.residues[f-1].atoms.size();l++){
                if(chain_p.residues[f-1].atoms[l].name==" C4'") f_1_C4=l;
            }
            RotVec[0]=chain_p.residues[f].atoms[f_C4].xyz[0]-chain_p.residues[f-1].atoms[f_1_C4].xyz[0];
            RotVec[1]=chain_p.residues[f].atoms[f_C4].xyz[1]-chain_p.residues[f-1].atoms[f_1_C4].xyz[1];
            RotVec[2]=chain_p.residues[f].atoms[f_C4].xyz[2]-chain_p.residues[f-1].atoms[f_1_C4].xyz[2];
        }else{
            vector<ResidueUnit>().swap(new_resi);
            return false;
        }

        for(int j=1;j<=sequence.size()-chain_p.residues[f].resi;j++){
            while(1){
                costheta=Gaussian(engine);
                if(costheta<=1&&costheta>=0) break;
            }
            phi=M_PI+M_PI*interval(engine);
            int f_1_C4=0;int f_2_C4=0;
            for(int l=0;l<new_resi[new_resi.size()-1].atoms.size();l++){
                if(new_resi[new_resi.size()-1].atoms[l].name==" C4'") f_1_C4=l;
            }
            new_coor_Rot(new_resi[new_resi.size()-1].atoms[f_1_C4].xyz,6.14,costheta,phi,atom.xyz,RotVec);
            resi_last.resi=chain_p.residues[f].resi+j;
            resi_last.resn[2]=toupper(sequence[resi_last.resi-1]);
            resi_last.atoms[0].xyz.swap(atom.xyz);
            new_resi.push_back(resi_last);
            for(int l=0;l<new_resi[new_resi.size()-1].atoms.size();l++){
                if(new_resi[new_resi.size()-1].atoms[l].name==" C4'") f_1_C4=l;
            }
            for(int l=0;l<new_resi[new_resi.size()-2].atoms.size();l++){
                if(new_resi[new_resi.size()-2].atoms[l].name==" C4'") f_2_C4=l;
            }
            RotVec[0]=new_resi[new_resi.size()-1].atoms[f_1_C4].xyz[0]-new_resi[new_resi.size()-2].atoms[f_2_C4].xyz[0];
            RotVec[1]=new_resi[new_resi.size()-1].atoms[f_1_C4].xyz[1]-new_resi[new_resi.size()-2].atoms[f_2_C4].xyz[1];
            RotVec[2]=new_resi[new_resi.size()-1].atoms[f_1_C4].xyz[2]-new_resi[new_resi.size()-2].atoms[f_2_C4].xyz[2];
        }

        vector<ResidueUnit> temp_resis;
        temp_resis.swap(new_resi);
        for(int j=1;j<temp_resis.size();j++){
            new_resi.push_back(temp_resis[j]);
        }
        vector<ResidueUnit>().swap(temp_resis);

        standardize_residue_order(new_resi);

        vector<double>().swap(RotVec);
        vector<AtomUnit>().swap(resi_new.atoms);
        vector<AtomUnit>().swap(resi_last.atoms);
        vector<double>().swap(atom.xyz);
        return true;
    }
    else if(f==-1&&r!=-1){
        resi_new.resi=chain_p.residues[r].resi;
        resi_new.resn=chain_p.residues[r].resn;
        int r_C4=0;
        for(int j=0;j<chain_p.residues[r].atoms.size();j++){
            if(chain_p.residues[r].atoms[j].name==" C4'") r_C4=j;
        }
        resi_new.atoms[0].xyz[0]=chain_p.residues[r].atoms[r_C4].xyz[0];
        resi_new.atoms[0].xyz[1]=chain_p.residues[r].atoms[r_C4].xyz[1];
        resi_new.atoms[0].xyz[2]=chain_p.residues[r].atoms[r_C4].xyz[2];
        new_resi.push_back(resi_new);

        if(r!=chain_p.residues.size()-1&&chain_p.residues[r+1].resi-chain_p.residues[r].resi==1){
            int r_1_C4=0;
            for(int l=0;l<chain_p.residues[r+1].atoms.size();l++){
                if(chain_p.residues[r+1].atoms[l].name==" C4'") r_1_C4=l;
            }
            RotVec[0]=chain_p.residues[r].atoms[r_C4].xyz[0]-chain_p.residues[r+1].atoms[r_1_C4].xyz[0];
            RotVec[1]=chain_p.residues[r].atoms[r_C4].xyz[1]-chain_p.residues[r+1].atoms[r_1_C4].xyz[1];
            RotVec[2]=chain_p.residues[r].atoms[r_C4].xyz[2]-chain_p.residues[r+1].atoms[r_1_C4].xyz[2];
        }else{
            vector<ResidueUnit>().swap(new_resi);
            return false;
        }


        for(int j=1;j<chain_p.residues[r].resi;j++){
            while(1){
                costheta=Gaussian(engine);
                if(costheta<=1&&costheta>=0) break;
            }
            phi=M_PI+M_PI*interval(engine);
            int r_1_C4=0;int r_2_C4=0;
            for(int l=0;l<new_resi[new_resi.size()-1].atoms.size();l++){
                if(new_resi[new_resi.size()-1].atoms[l].name==" C4'") r_1_C4=l;
            }
            //cout<<new_resi[new_resi.size()-1].resi<<"\t"
            //    <<new_resi[new_resi.size()-1].atoms[r_1_C4].xyz[0]<<"\t"
            //    <<new_resi[new_resi.size()-1].atoms[r_1_C4].xyz[1]<<"\t"
            //    <<new_resi[new_resi.size()-1].atoms[r_1_C4].xyz[2]<<endl;
            cout<<endl;
            new_coor_Rot(new_resi[new_resi.size()-1].atoms[r_1_C4].xyz,6.14,costheta,phi,atom.xyz,RotVec);
            cout<<atom.xyz[0]<<"\t"
                <<atom.xyz[1]<<"\t"
                <<atom.xyz[2]<<endl;
            resi_last.resi=chain_p.residues[r].resi-j;
            resi_last.resn[2]=toupper(sequence[resi_last.resi-1]);
            resi_last.atoms[0].xyz.swap(atom.xyz);
            new_resi.push_back(resi_last);
            for(int l=0;l<new_resi[new_resi.size()-1].atoms.size();l++){
                if(new_resi[new_resi.size()-1].atoms[l].name==" C4'") r_1_C4=l;
            }
            for(int l=0;l<new_resi[new_resi.size()-2].atoms.size();l++){
                if(new_resi[new_resi.size()-2].atoms[l].name==" C4'") r_2_C4=l;
            }
            RotVec[0]=new_resi[new_resi.size()-1].atoms[r_1_C4].xyz[0]-new_resi[new_resi.size()-2].atoms[r_2_C4].xyz[0];
            RotVec[1]=new_resi[new_resi.size()-1].atoms[r_1_C4].xyz[1]-new_resi[new_resi.size()-2].atoms[r_2_C4].xyz[1];
            RotVec[2]=new_resi[new_resi.size()-1].atoms[r_1_C4].xyz[2]-new_resi[new_resi.size()-2].atoms[r_2_C4].xyz[2];
        }
        vector<ResidueUnit> temp_resis;
        temp_resis.swap(new_resi);
        for(int j=1;j<temp_resis.size();j++){
            //cout<<temp_resis[j].atoms[0].xyz[0]<<"\t"<<temp_resis[j].atoms[0].xyz[1]<<"\t"<<temp_resis[j].atoms[0].xyz[2]<<endl;
            new_resi.push_back(temp_resis[j]);
        }
        vector<ResidueUnit>().swap(temp_resis);

        standardize_residue_order(new_resi);

        vector<double>().swap(RotVec);
        vector<AtomUnit>().swap(resi_new.atoms);
        vector<AtomUnit>().swap(resi_last.atoms);
        vector<double>().swap(atom.xyz);
        return true;
    }
    else{
        vector<double>().swap(RotVec);
        vector<AtomUnit>().swap(resi_new.atoms);
        vector<AtomUnit>().swap(resi_last.atoms);
        vector<double>().swap(atom.xyz);
        return false;
    }
}

bool fill_SRS(string sequence, ChainUnit &chain_p,int f,int r,vector<ResidueUnit> &new_resi){
    if(chain_p.residues[r].resi<chain_p.residues[f].resi){
        cout<<"invalid input!"<<endl;
        return false;
    }
    vector<double> RotVec1;
    vector<double> RotVec2;
    RotVec1.assign(3,0);
    RotVec2.assign(3,0);
    RotVec1[2]=1;
    RotVec2[2]=1;
    

    ResidueUnit resi_n;
    ResidueUnit resi_nf;
    ResidueUnit resi_nr;
    resi_n.het=false;
    resi_n.icode=' ';
    resi_n.resn="   ";
    resi_nf.het=false;
    resi_nf.icode=' ';
    resi_nf.resn="   ";
    resi_nr.het=false;
    resi_nr.icode=' ';
    resi_nr.resn="   ";
    AtomUnit atom;
    atom.movable=1;
    atom.name=chain_p.residues[f].atoms[5].name;
    atom.xyz.assign(3,0);
    resi_n.atoms.assign(1,atom);
    resi_nf.atoms.assign(1,atom);
    resi_nr.atoms.assign(1,atom);
    normal_distribution<double> Gaussian(0.848,0.123);
    uniform_real_distribution<double> interval(-1.0,1.0);
    default_random_engine engine;
    engine.seed(4922995424);
    //engine.seed(time(NULL));
    //cout<<"SRS"<<endl;
    double costheta;
    double phi;
    vector<ResidueUnit> res_f;
    vector<ResidueUnit> res_r;
    resi_nf.resi=chain_p.residues[f].resi;
    resi_nf.resn=chain_p.residues[f].resn;
    int f_C4=0;
    for(int i=0;i<chain_p.residues[f].atoms.size();i++){
        if(chain_p.residues[f].atoms[i].name==" C4'") f_C4=i;
    } //There must be a C4' in the residue
    resi_nf.atoms[0].xyz[0]=chain_p.residues[f].atoms[f_C4].xyz[0];
    resi_nf.atoms[0].xyz[1]=chain_p.residues[f].atoms[f_C4].xyz[1];
    resi_nf.atoms[0].xyz[2]=chain_p.residues[f].atoms[f_C4].xyz[2];
    res_f.push_back(resi_nf);
    resi_nr.resi=chain_p.residues[r].resi;
    resi_nr.resn=chain_p.residues[r].resn;
    int r_C4=0;
    for(int i=0;i<chain_p.residues[r].atoms.size();i++){
        if(chain_p.residues[r].atoms[i].name==" C4'") r_C4=i;
    }
    resi_nr.atoms[0].xyz[0]=chain_p.residues[r].atoms[r_C4].xyz[0];
    resi_nr.atoms[0].xyz[1]=chain_p.residues[r].atoms[r_C4].xyz[1];
    resi_nr.atoms[0].xyz[2]=chain_p.residues[r].atoms[r_C4].xyz[2];
    res_r.push_back(resi_nr);
    int p=40;
    cout<<Points2Distance(resi_nf.atoms[0].xyz,resi_nr.atoms[0].xyz)<<endl;

    while(p)
    {
    //cout<<"SRS"<<endl;
    int t=0;

    if(f!=0&&chain_p.residues[f].resi-chain_p.residues[f-1].resi==1){
        int f_1_C4=0;
        for(int j=0;j<chain_p.residues[f-1].atoms.size();j++){
            if(chain_p.residues[f-1].atoms[j].name==" C4'") f_1_C4=j;
        }
        RotVec1[0]=chain_p.residues[f].atoms[f_C4].xyz[0]-chain_p.residues[f-1].atoms[f_1_C4].xyz[0];
        RotVec1[1]=chain_p.residues[f].atoms[f_C4].xyz[1]-chain_p.residues[f-1].atoms[f_1_C4].xyz[1];
        RotVec1[2]=chain_p.residues[f].atoms[f_C4].xyz[2]-chain_p.residues[f-1].atoms[f_1_C4].xyz[2];
        //cout<<"i1"<<"\t"<<RotVec1[0]<<"\t"<<RotVec1[1]<<"\t"<<RotVec1[2]<<endl;
    }
    if(r!=chain_p.residues.size()-1&&chain_p.residues[r+1].resi-chain_p.residues[r].resi==1){
        int r_1_C4=0;
        for(int j=0;j<chain_p.residues[r+1].atoms.size();j++){
            if(chain_p.residues[r+1].atoms[j].name==" C4'") r_1_C4=j;
        }
        RotVec2[0]=chain_p.residues[r].atoms[r_C4].xyz[0]-chain_p.residues[r+1].atoms[r_1_C4].xyz[0];
        RotVec2[1]=chain_p.residues[r].atoms[r_C4].xyz[1]-chain_p.residues[r+1].atoms[r_1_C4].xyz[1];
        RotVec2[2]=chain_p.residues[r].atoms[r_C4].xyz[2]-chain_p.residues[r+1].atoms[r_1_C4].xyz[2];
        //cout<<"i2"<<"\t"<<RotVec2[0]<<"\t"<<RotVec2[1]<<"\t"<<RotVec2[2]<<endl;
    }

    for(int i=1;i<chain_p.residues[r].resi-chain_p.residues[f].resi;i++){
        double p_atom=false;
        while(!p_atom&&t<20*(chain_p.residues[r].resi-chain_p.residues[f].resi)){
            //cout<<i<<"\t"<<res_f.size()<<"\t"<<res_r.size()<<"\t"<<t<<endl;
        if(i%2){
            while(1){
                costheta=Gaussian(engine);
                if(costheta<=1&&costheta>=0) break;
            }
            //costheta=interval(engine);//;inter()
            phi=M_PI+M_PI*interval(engine);//inter();
            int f_1_C4=0;
            int f_2_C4=0;
            for(int k=0;k<res_f[res_f.size()-1].atoms.size();k++){
                if(res_f[res_f.size()-1].atoms[k].name==" C4'") f_1_C4=k;
            }
            new_coor_Rot(res_f[res_f.size()-1].atoms[f_1_C4].xyz,6.14,costheta,phi,atom.xyz,RotVec1);
            t++;
            if(ok_SRS(res_f,res_r,atom.xyz,i)){
                //cout<<"costheta+phi"<<"\t"<<costheta<<"\t"<<phi<<endl;
                //cout<<0.5<<"\t"<<RotVec1[0]<<"\t"<<RotVec1[1]<<"\t"<<RotVec1[2]<<endl;```````````
                p_atom=true;
                resi_n.resi=chain_p.residues[f].resi+i/2+1;
                resi_n.resn[2]=toupper(sequence[resi_n.resi-1]);
                resi_n.atoms[0].xyz.swap(atom.xyz);
                res_f.push_back(resi_n);
                for(int l=0;l<res_f[res_f.size()-1].atoms.size();l++){
                    if(res_f[res_f.size()-1].atoms[l].name==" C4'") f_1_C4=l;
                }
                for(int l=0;l<res_f[res_f.size()-1].atoms.size();l++){
                    if(res_f[res_f.size()-2].atoms[l].name==" C4'") f_2_C4=l;
                }
                RotVec1[0]=res_f[res_f.size()-1].atoms[f_1_C4].xyz[0]-res_f[res_f.size()-2].atoms[f_2_C4].xyz[0];
                RotVec1[1]=res_f[res_f.size()-1].atoms[f_1_C4].xyz[1]-res_f[res_f.size()-2].atoms[f_2_C4].xyz[1];
                RotVec1[2]=res_f[res_f.size()-1].atoms[f_1_C4].xyz[2]-res_f[res_f.size()-2].atoms[f_2_C4].xyz[2];
                //cout<<1<<"\t"<<RotVec1[0]<<"\t"<<RotVec1[1]<<"\t"<<RotVec1[2]<<endl;
            }
            
        }
        else{
            while(1){
                costheta=Gaussian(engine);
                if(costheta<=1&&costheta>=0) break;
            }
            //costheta=interval(engine);//inter();
            phi=M_PI+M_PI*interval(engine);//inter();
            int r_1_C4=0;
            int r_2_C4=0;
            for(int k=0;k<res_r[res_r.size()-1].atoms.size();k++){
                if(res_r[res_r.size()-1].atoms[k].name==" C4'") r_1_C4=k;
            }
            new_coor_Rot(res_r[res_r.size()-1].atoms[r_1_C4].xyz,6.14,costheta,phi,atom.xyz,RotVec2);
            t++;
            if(ok_SRS(res_f,res_r,atom.xyz,i)){
                //cout<<"costheta+phi"<<"\t"<<costheta<<"\t"<<phi<<endl;
                //cout<<1.5<<"\t"<<RotVec2[0]<<"\t"<<RotVec2[1]<<"\t"<<RotVec2[2]<<endl;
                p_atom=true;
                resi_n.resi=chain_p.residues[r].resi-i/2;
                resi_n.resn[2]=toupper(sequence[resi_n.resi-1]);
                resi_n.atoms[0].xyz.swap(atom.xyz);
                res_r.push_back(resi_n);
                for(int l=0;l<res_r[res_r.size()-1].atoms.size();l++){
                    if(res_r[res_r.size()-1].atoms[l].name==" C4'") r_1_C4=l;
                }
                for(int l=0;l<res_r[res_r.size()-1].atoms.size();l++){
                    if(res_r[res_r.size()-2].atoms[l].name==" C4'") r_2_C4=l;
                }
                RotVec2[0]=res_r[res_r.size()-1].atoms[r_1_C4].xyz[0]-res_r[res_r.size()-2].atoms[r_2_C4].xyz[0];
                RotVec2[1]=res_r[res_r.size()-1].atoms[r_1_C4].xyz[1]-res_r[res_r.size()-2].atoms[r_2_C4].xyz[1];
                RotVec2[2]=res_r[res_r.size()-1].atoms[r_1_C4].xyz[2]-res_r[res_r.size()-2].atoms[r_2_C4].xyz[2];
                //cout<<2<<"\t"<<RotVec2[0]<<"\t"<<RotVec2[1]<<"\t"<<RotVec2[2]<<endl;
            }
            
        }
        }
        if(t<20*(chain_p.residues[r].resi-chain_p.residues[f].resi)) cout<<p<<"\t"<<i<<endl;
    }
    //cout<<t<<"\t"<<resi_r.resi-resi_f.resi<<endl;
    int f_1_C4=0;
    int r_1_C4=0;
    for(int i=0;i<res_f[res_f.size()-1].atoms.size();i++){
        if(res_f[res_f.size()-1].atoms[i].name==" C4'") f_1_C4=i;
    }
    for(int i=0;i<res_r[res_r.size()-1].atoms.size();i++){
        if(res_r[res_r.size()-1].atoms[i].name==" C4'") r_1_C4=i;
    }
    if(t<20*(chain_p.residues[r].resi-chain_p.residues[f].resi)) p=0;
    else if(Points2Distance(res_f[res_f.size()-1].atoms[f_1_C4].xyz,res_r[res_r.size()-1].atoms[r_1_C4].xyz)>5.6*(res_r[res_r.size()-1].resi-res_f[res_f.size()-1].resi)
          &&Points2Distance(res_f[res_f.size()-1].atoms[f_1_C4].xyz,res_r[res_r.size()-1].atoms[r_1_C4].xyz)<6.0*(res_r[res_r.size()-1].resi-res_f[res_f.size()-1].resi)){
        p=0;
        int gap_atom_num=res_r[res_r.size()-1].resi-res_f[res_f.size()-1].resi-1;
        double x_distance=res_r[res_r.size()-1].atoms[r_1_C4].xyz[0]-res_f[res_f.size()-1].atoms[f_1_C4].xyz[0];
        double y_distance=res_r[res_r.size()-1].atoms[r_1_C4].xyz[1]-res_f[res_f.size()-1].atoms[f_1_C4].xyz[1];
        double z_distance=res_r[res_r.size()-1].atoms[r_1_C4].xyz[2]-res_f[res_f.size()-1].atoms[f_1_C4].xyz[2];
        for(int j=0;j<gap_atom_num;j++){
            resi_n.resi=res_f[res_f.size()-1].resi+1;
            resi_n.resn[2]=toupper(sequence[resi_n.resi-1]);
            resi_n.atoms[0].xyz[0]=res_f[res_f.size()-1].atoms[f_1_C4].xyz[0]+x_distance/(gap_atom_num+1);
            resi_n.atoms[0].xyz[1]=res_f[res_f.size()-1].atoms[f_1_C4].xyz[1]+y_distance/(gap_atom_num+1);
            resi_n.atoms[0].xyz[2]=res_f[res_f.size()-1].atoms[f_1_C4].xyz[2]+z_distance/(gap_atom_num+1);
            res_f.push_back(resi_n);
        }
    }
    else{
        vector<ResidueUnit>().swap(res_f);
        vector<ResidueUnit>().swap(res_r);
        res_f.push_back(resi_nf);
        res_r.push_back(resi_nr);
        p--;
        if(p==0) return false;
    }
    }
    for(int i=1;i<res_f.size();i++) new_resi.push_back(res_f[i]);
    for(int i=1;i<res_r.size();i++) new_resi.push_back(res_r[res_r.size()-i]);
    vector<ResidueUnit>().swap(res_f);
    vector<ResidueUnit>().swap(res_r);
    
    return true;
}

void ChainRotation(vector<ResidueUnit> &chain, vector<double> &axisA, vector<double> &axisB, double angle){
    vector<vector<double>> coordinates;
    for(size_t i=0;i<chain.size();i++){
        for(int j=0;j<chain[i].atoms.size();j++){
            coordinates.push_back(chain[i].atoms[j].xyz);
        }
    }
    GroupRotation(axisA,axisB,angle,coordinates);
    int k=0;
    for(size_t i=0;i<chain.size();i++){
        for(int j=0;j<chain[i].atoms.size();j++){
            chain[i].atoms[j].xyz.swap(coordinates[k]);
            k++;
        }
    }
}

bool SingleRotation(const vector<double> &axisA, const vector<double> &axisB, double angle, vector<double> &point)
{
   if(axisA.size()!=3 || axisB.size()!=3)
   {
      cerr<<"Error in SingleRotation()"<<endl;
      return false;
   }
   
   vector<double> axis(3, 0);
   axis[0]=axisB[0]-axisA[0]; 
   axis[1]=axisB[1]-axisA[1];
   axis[2]=axisB[2]-axisA[2]; 
   
   matrix rotmtx;
   if(!RotationMatrixB(axis, deg2rad(angle), rotmtx)) return false;
   
   //vector<double> point_A(4, 0);
   double point_A[3];

      point_A[0]=point[0]-axisA[0];
      point_A[1]=point[1]-axisA[1];
      point_A[2]=point[2]-axisA[2];

      point=axisA;
      point[0]+=rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2]+rotmtx[0][3];
      point[1]+=rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2]+rotmtx[1][3];
      point[2]+=rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2]+rotmtx[2][3];
   return true;
}

/*vector<vector<ResidueUnit> > new_loop_rotation(vector<ResidueUnit> &new_resi, vector<double> xyz_f,vector<double> xyz_r, int times){
    vector<vector<ResidueUnit> > multi_loops;
    multi_loops.push_back(new_resi);
    for(int i=0;i<times;i++){
        int a1=JKISS()%(new_resi.size()+2);
        int a2=JKISS()%(new_resi.size()+2);
        cout<<a1<<" "<<a2<<endl;
        while(abs(a2-a1)<2) a2=JKISS()%(new_resi.size()+2);
        if(a1>a2){
            int a;
            a=a1;a1=a2;a2=a;
        }
        for(int j=a1;j<a2-1;j++){
            for(int k=0;k<new_resi[j].atoms.size();k++){
                if(a1==0&&a2==new_resi.size()+1){
                    SingleRotation(xyz_f,xyz_r,10,new_resi[j].atoms[k].xyz);
                }else if(a1==0){
                    SingleRotation(xyz_f,new_resi[a2-1].atoms[0].xyz,10,new_resi[j].atoms[k].xyz);
                }else if(a2==new_resi.size()+1){
                    SingleRotation(new_resi[a1-1].atoms[0].xyz,xyz_r,10,new_resi[j].atoms[k].xyz);
                }else{
                    SingleRotation(new_resi[a1-1].atoms[0].xyz,new_resi[a2-1].atoms[0].xyz,10,new_resi[j].atoms[k].xyz);
                }
            }
        }
        multi_loops.push_back(new_resi);
    }
    return multi_loops;
}

vector<int> atom_clashes_C4_quick(ModelUnit &pep, vector<vector<ResidueUnit> > &new_resis, ResidueUnit resi_f, ResidueUnit resi_r){
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
    int clash;
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
    return clashes;
}*/


double loop_value(double* model_array,size_t model_size,double* loop_array,int loop_size){
    double value=0;
    double x,y,z;
    
    for(size_t i=0;i<model_size;i++){
        for(int j=1;j<=loop_size;j++){
            x=model_array[3*i+0]-loop_array[3*j+0];
            y=model_array[3*i+1]-loop_array[3*j+1];
            z=model_array[3*i+2]-loop_array[3*j+2];
            if(x*x+y*y+z*z<38.44) value+=(38.44-x*x-y*y-z*z);
        }
    }
    return value;
}

bool rotate_loop(double* loop_array,double* temp_array,int loop_size){
	int a=JKISS()%(loop_size+2);
	int b=JKISS()%(loop_size+2);
	while(abs(b-a)<2) b=JKISS()%(loop_size+2);
	int head,tail;
	if(a<b){
		head=a;
		tail=b;
	}else{
		head=b;
		tail=a;
	}
    
    vector<double> axis(3, 0);
    axis[0]=loop_array[3*b+0]-loop_array[3*a+0]; 
    axis[1]=loop_array[3*b+1]-loop_array[3*a+1]; 
    axis[2]=loop_array[3*b+2]-loop_array[3*a+2];

    matrix rotmtx;
   if(!RotationMatrixB(axis, deg2rad(12), rotmtx)) return false; 

   double point[3];
   for(int i=a+1; i<b; i++)
   {
      point[0]=loop_array[3*i+0]-loop_array[3*a+0];
      point[1]=loop_array[3*i+1]-loop_array[3*a+1];
      point[2]=loop_array[3*i+2]-loop_array[3*a+2];

      temp_array[3*i+0]=loop_array[3*a+0];
      temp_array[3*i+1]=loop_array[3*a+1];
      temp_array[3*i+2]=loop_array[3*a+2];
      temp_array[3*i+0]+=rotmtx[0][0]*point[0]+rotmtx[0][1]*point[1]+rotmtx[0][2]*point[2]+rotmtx[0][3];
      temp_array[3*i+1]+=rotmtx[1][0]*point[0]+rotmtx[1][1]*point[1]+rotmtx[1][2]*point[2]+rotmtx[1][3];
      temp_array[3*i+2]+=rotmtx[2][0]*point[0]+rotmtx[2][1]*point[1]+rotmtx[2][2]*point[2]+rotmtx[2][3];
   }
   return true;
}

inline bool second_is_better(double first, double second){
    if(second<first) return true;
    else return false;
}

bool select_loop(ModelUnit &model, vector<ResidueUnit> &new_resi, vector<double> &axis_f, vector<double> &axis_r){
	
	//prepare loop
    int loop_size=new_resi.size();
    if(loop_size==0) return false;
	int loop_start=new_resi[0].resi;
	
	double* loop_array=new double[3*(loop_size+2)];
	for(int i=1;i<=loop_size;i++){
		loop_array[3*i+0]=new_resi[i-1].atoms[0].xyz[0];
		loop_array[3*i+1]=new_resi[i-1].atoms[0].xyz[1];
		loop_array[3*i+2]=new_resi[i-1].atoms[0].xyz[2];
	}
	loop_array[0]=axis_f[0];
	loop_array[1]=axis_f[1];
	loop_array[2]=axis_f[2];
	loop_array[3*(loop_size+1)+0]=axis_r[0];
	loop_array[3*(loop_size+1)+1]=axis_r[1];
	loop_array[3*(loop_size+1)+2]=axis_r[2];
	
	//prepare model
	vector<vector<double> > atoms_concerned;

    for(int i=0;i<model.chains.size();i++){
        for(size_t j=0;j<model.chains[i].residues.size();j++){
            for(int k=0;k<model.chains[i].residues[j].atoms.size();k++){
                if(Points2Distance(axis_f,model.chains[i].residues[j].atoms[k].xyz)
                +Points2Distance(axis_r,model.chains[i].residues[j].atoms[k].xyz)
                <6.2*(loop_size+1)+2*6.2)
                atoms_concerned.push_back(model.chains[i].residues[j].atoms[k].xyz);
            }
        }
    }


	size_t model_size=atoms_concerned.size();
	double* model_array=new double[3*model_size];
	for(int i=0;i<model_size;i++){
		model_array[3*i+0]=atoms_concerned[i][0];
		model_array[3*i+1]=atoms_concerned[i][1];
		model_array[3*i+2]=atoms_concerned[i][2];
	}
	
	//rotate
	double E_value=loop_value(model_array,model_size,loop_array,loop_size);
	double temp_value;
	double* temp_array=new double[3*(loop_size+2)];
    cout<<model_size<<endl;
    cout<<E_value<<endl;
    cout<<loop_size<<endl;
    for(int i=0;i<3*(loop_size+2);i++) temp_array[i]=loop_array[i];
	for(int i=0;i<200;i++){
		//rotate loop
		rotate_loop(loop_array,temp_array,loop_size);
        
		//calculate value
		temp_value=loop_value(model_array,model_size,temp_array,loop_size);
		//compare and adjust
		if(second_is_better(E_value,temp_value)){
            E_value=temp_value;
            for(int j=0;j<3*(loop_size+2);j++) loop_array[j]=temp_array[j];
        }else for(int j=0;j<3*(loop_size+2);j++) temp_array[j]=loop_array[j];
	}
	for(int i=1;i<=loop_size;i++){
		new_resi[i-1].atoms[0].xyz[0]=loop_array[3*i+0];
		new_resi[i-1].atoms[0].xyz[1]=loop_array[3*i+1];
		new_resi[i-1].atoms[0].xyz[2]=loop_array[3*i+2];
	}
    delete[] model_array;
    delete[] loop_array;
    delete[] temp_array;
    return true;
}


int multi_loops_SRS(const string sequence, ModelUnit &model_p, ChainUnit &chain_p){
    vector<vector<ResidueUnit> > gapif;
    sizet_sizet_vc gaps=gapinfo(sequence,chain_p,gapif,0,1);
    vector<ResidueUnit> new_resi;
    vector<vector<ResidueUnit> > multi_loops;
    sizet_sizet_vc gaps_temp;
    cout<<gaps.size()<<endl;
    for(size_t g=0;g<gaps.size();g++){
        int f,r;
        f=-1;r=-1;
        for(size_t j=0;j<chain_p.residues.size();j++){
            if(chain_p.residues[j].resi==gaps[g].first-1) f=j;
            if(chain_p.residues[j].resi==gaps[g].second) r=j;
        }
        /*for(int k=0;k<5;k++){
            if(Points2Distance(chain_p.residues[f].atoms[5].xyz,chain_p.residues[r].atoms[5].xyz)>5.9*(chain_p.residues[r].resi-chain_p.residues[f].resi)) fill_linear(sequence,chain_p.residues[f],chain_p.residues[r],new_resi);
            else fill_SRS(sequence,chain_p.residues[f],chain_p.residues[r],new_resi);
            multi_loops.push_back(new_resi);
            vector<ResidueUnit>().swap(new_resi);
        }*/
        cout<<gaps[g].first-1<<"\t"<<gaps[g].second<<endl;
        cout<<f<<"\t"<<r<<endl;
        if(f==-1||r==-1){
            fill_SRS_terminal(sequence,chain_p,f,r,new_resi);
            cout<<new_resi.size()<<endl;
            fill_new_residue(chain_p,new_resi);
            vector<ResidueUnit>().swap(new_resi);
        } 
        else if(Points2Distance(chain_p.residues[f].atoms[5].xyz,chain_p.residues[r].atoms[5].xyz)>5.6*(chain_p.residues[r].resi-chain_p.residues[f].resi)||chain_p.residues[r].resi-chain_p.residues[f].resi<3){
            fill_linear(sequence,chain_p.residues[f],chain_p.residues[r],new_resi);
            fill_new_residue(chain_p,new_resi);
            vector<ResidueUnit>().swap(new_resi);
        } 
        else{
            cout<<"AAA"<<endl;
            if(fill_SRS(sequence,chain_p,f,r,new_resi)){
                select_loop(model_p,new_resi,chain_p.residues[f].atoms[5].xyz,chain_p.residues[r].atoms[5].xyz);
                fill_new_residue(chain_p,new_resi);
                vector<ResidueUnit>().swap(new_resi);
            }
            else{
                cout<<"SRS failed, using linear fill instead."<<endl;
                fill_linear(sequence,chain_p.residues[f],chain_p.residues[r],new_resi);
                fill_new_residue(chain_p,new_resi);
                vector<ResidueUnit>().swap(new_resi);
            }
            //cout<<"SRS"<<new_resi.size()<<endl;

            /*vector<vector<ResidueUnit> > new_resis=new_loop_rotation(new_resi,chain_p.residues[f].atoms[5].xyz,chain_p.residues[r].atoms[5].xyz,15);
            cout<<"loop"<<endl;
            vector<ResidueUnit> best_resi=atom_clashes_C4_quick(model_p,new_resis,chain_p.residues[f],chain_p.residues[r]);
            cout<<"clash"<<endl;*/
            //vector<vector<ResidueUnit> >().swap(new_resis);
            //vector<ResidueUnit>().swap(best_resi);
        } 
        cout<<"find clashes"<<endl;
        //cout<<"out"<<endl;
        
        standardize_residue_order(chain_p);
        cout<<chain_p.residues[0].resi<<endl;
    }
    return 1;
}

int fill_SRS(const string sequence,ChainUnit &chain_p){
    sizet_sizet_vc gaps;
    size_t start=1;
    for(size_t i=0;i<chain_p.residues.size();i++){
        if(chain_p.residues[i].resi==start){
            start++;
        }else{
            gaps.push_back(pair<size_t,size_t>(start,chain_p.residues[i].resi));
            start=chain_p.residues[i].resi+1;
        }
    }
    if(start<sequence.size()+1) gaps.push_back(pair<size_t,size_t>(start,sequence.size()+1));
    vector<ResidueUnit> new_resi;
    sizet_sizet_vc gaps_temp;
    for(size_t g=0;g<gaps.size();g++){
        bool finish;
        size_t f,r;
        for(size_t j=0;j<chain_p.residues.size();j++){
            if(chain_p.residues[j].resi==gaps[g].first-1) f=j;
            if(chain_p.residues[j].resi==gaps[g].second) r=j;
        }
        finish=fill_SRS(sequence,chain_p,f,r,new_resi);
        //cout<<"out"<<endl;
        for(size_t j=0;j<new_resi.size();j++){
            chain_p.residues.push_back(new_resi[j]);
        }
        standardize_residue_order(chain_p);
        vector<ResidueUnit>().swap(new_resi);
        //cout<<"out"<<endl;
    }
    return 1;
}

/*FLORA_7 - Spherical Random Sampling*/
void fill_SRS(Model_F &model_f,ModelUnit &model_p){
    int i1;
    for(size_t i=0;i<model_p.chains.size();i++){
        //i1=fill_SRS(model_f.chains_F[i],model_p.chains[i]);
        i1=multi_loops_SRS(model_f.chains_F[i],model_p,model_p.chains[i]);
    }
    standardize_residue_order(model_p);
}

/*FLORA_0 - linear filling*/
void fill_gap_whole(Model_F &model_f,ModelUnit &model_p, const int option=0){
    int i1;
    for(size_t i=0;i<model_f.chains_F.size();i++){
        if(option==0) i1=fill_linear(model_f.chains_F[i],model_p.chains[i]);
    }
    standardize_residue_order(model_p);
}

bool RMSD_raw(pair<size_t,size_t> gap_if,ChainUnit &chain_in,ChainUnit &chain_ref,int option){
    ofstream fp_out;
    fp_out.open("RMSD-r.txt",ios::app);
    int gap_length=gap_if.second-gap_if.first;
    int front=-7;
    int rear=-7;
    double gap_distance=0;
    for(int i=0;i<chain_in.residues.size();i++){
        if(chain_in.residues[i].resi==gap_if.first&&i>0) front=i-1;
        if(chain_in.residues[i].resi==gap_if.second) rear=i; 
    }
    if(front!=-7&&rear!=-7) gap_distance=Points2Distance(chain_in.residues[front].atoms[5].xyz,chain_in.residues[rear].atoms[5].xyz);
    vector<ResidueUnit> gap_in;
    vector<ResidueUnit> gap_ref;
    for(int i=0;i<chain_in.residues.size();i++){
        if(chain_in.residues[i].resi==gap_if.first){
            for(int j=0;j<gap_length;j++){
                gap_in.push_back(chain_in.residues[i+j]);
            }
        }
    }
    for(int i=0;i<chain_ref.residues.size();i++){
        if(chain_ref.residues[i].resi==gap_if.first){
            for(int j=0;j<gap_length;j++){
                gap_ref.push_back(chain_ref.residues[i+j]);
            }
        }
    }
    if(gap_in.size()!=gap_ref.size()){
        cout<<"Residues Unmatched!"<<endl;
        fp_out.close();
        return false;
    }
    double AT_rate=0;
    for(int i=0;i<gap_in.size();i++){
        if(gap_in[i].resn=="  A"||gap_in[i].resn=="  U") AT_rate+=1;
    }
    AT_rate/=gap_in.size();

    double count=0;
    double sd=0;
    for(int i=0;i<gap_in.size();i++){
        if(gap_in[i].atoms.size()!=gap_ref[i].atoms.size()){
            cout<<"Atoms Unmatched!"<<endl;
            fp_out.close();
            return false;
        } 
        for(int j=0;j<gap_in[i].atoms.size();j++){
            sd+=Points2Distance2(gap_ref[i].atoms[j].xyz,gap_in[i].atoms[j].xyz);
            count+=1;
        }
    }
    double MSD_value=sd/count;
    fp_out<<gap_length<<"\t"<<gap_distance<<"\t"<<AT_rate<<"\t"<<pow(MSD_value,0.5)<<"\t"<<option<<endl;

    fp_out.close();
    return true;
}

bool RMSD(pair<size_t,size_t> gap_if,ChainUnit &chain_in,ChainUnit &chain_ref,int option){
    ofstream fp_out;
    fp_out.open("RMSD-o.txt",ios::app);
    int gap_length=gap_if.second-gap_if.first;
    int front=-7;
    int rear=-7;
    double gap_distance=0;
    for(int i=0;i<chain_in.residues.size();i++){
        if(chain_in.residues[i].resi==gap_if.first&&i>0) front=i-1;
        if(chain_in.residues[i].resi==gap_if.second) rear=i; 
    }
    if(front!=-7&&rear!=-7) gap_distance=Points2Distance(chain_in.residues[front].atoms[5].xyz,chain_in.residues[rear].atoms[5].xyz);
    vector<ResidueUnit> gap_in;
    vector<ResidueUnit> gap_ref;
    for(int i=0;i<chain_in.residues.size();i++){
        if(chain_in.residues[i].resi==gap_if.first){
            for(int j=0;j<gap_length;j++){
                gap_in.push_back(chain_in.residues[i+j]);
            }
        }
    }
    for(int i=0;i<chain_ref.residues.size();i++){
        if(chain_ref.residues[i].resi==gap_if.first){
            for(int j=0;j<gap_length;j++){
                gap_ref.push_back(chain_ref.residues[i+j]);
            }
        }
    }
    if(gap_in.size()!=gap_ref.size()){
        cout<<"Residues Unmatched!"<<endl;
        fp_out.close();
        return false;
    }
    double AT_rate=0;
    for(int i=0;i<gap_in.size();i++){
        if(gap_in[i].resn=="  A"||gap_in[i].resn=="  U") AT_rate+=1;
    }
    AT_rate/=gap_in.size();

    double count=0;
    double sd=0;
    for(int i=0;i<gap_in.size();i++){
        if(gap_in[i].atoms.size()!=gap_ref[i].atoms.size()){
            cout<<"Atoms Unmatched!"<<endl;
            fp_out.close();
            return false;
        } 
        for(int j=0;j<gap_in[i].atoms.size();j++){
            sd+=Points2Distance2(gap_ref[i].atoms[j].xyz,gap_in[i].atoms[j].xyz);
            count+=1;
        }
    }
    double MSD_value=sd/count;
    fp_out<<gap_length<<"\t"<<gap_distance<<"\t"<<AT_rate<<"\t"<<pow(MSD_value,0.5)<<"\t"<<option<<endl;

    fp_out.close();
    return true;
}

bool can_handle(Model_F &fasta_in,ModelUnit &pdb_in,const int option=0){
    for(int i=0;i<fasta_in.chains_F.size();i++){
        sizet_sizet_vc gapif=gapinfo(fasta_in.chains_F[i],pdb_in.chains[i]);
        sizet_sizet_vc frends;
        for(int j=0;j<gapif.size();j++){
            pair<size_t,size_t> frend;
            if(j==0) frend.first=gapif[j].first-1;
            else frend.first=gapif[j].first-gapif[j-1].second;
            if(j==gapif.size()-1) frend.second=fasta_in.chains_F[i].size()+1-gapif[j].second;
            else frend.second=gapif[j+1].first-gapif[j].second;
            frends.push_back(frend);
        }
        if(option==0){
            for(int k=0;k<frends.size();k++){
                if(frends[k].first>=1&&frends[k].second>=1){
                    sizet_sizet_vc().swap(frends);
                    sizet_sizet_vc().swap(gapif);
                    return true;
                }
            }
        }
        else if(option==1){
            for(int k=0;k<frends.size();k++){
                if(frends[k].first>=4||frends[k].second>=4){
                    sizet_sizet_vc().swap(frends);
                    sizet_sizet_vc().swap(gapif);
                    return true;
                }
            }
        }
        else if(option==2){
            for(int k=0;k<frends.size();k++){
                if(frends[k].first+frends[k].second>=2){
                    sizet_sizet_vc().swap(frends);
                    sizet_sizet_vc().swap(gapif);
                    return true;
                }
            }
        }
        else if(option==3){
            for(int k=0;k<frends.size();k++){
                if(frends[k].first+frends[k].second>=2){
                    sizet_sizet_vc().swap(frends);
                    sizet_sizet_vc().swap(gapif);
                    return true;
                }
            }
        }
        else if(option==4){
            for(int k=0;k<frends.size();k++){
                if(frends[k].first+frends[k].second>=2){
                    sizet_sizet_vc().swap(frends);
                    sizet_sizet_vc().swap(gapif);
                    return true;
                }
            }
        }
        else if(option==5){
            for(int k=0;k<frends.size();k++){
                if(frends[k].first+frends[k].second>=2){
                    sizet_sizet_vc().swap(frends);
                    sizet_sizet_vc().swap(gapif);
                    return true;
                }
            }
        }
        else if(option==6){
            for(int k=0;k<frends.size();k++){
                if(frends[k].first>=1&&frends[k].second>=1){
                    sizet_sizet_vc().swap(frends);
                    sizet_sizet_vc().swap(gapif);
                    return true;
                }
            }
        }
        else if(option==7){
            for(int k=0;k<frends.size();k++){
                if(frends[k].first+frends[k].second>=2){
                    sizet_sizet_vc().swap(frends);
                    sizet_sizet_vc().swap(gapif);
                    return true;
                }
            }
        }
        sizet_sizet_vc().swap(frends);
        sizet_sizet_vc().swap(gapif);
    }
    return false;
}

/*choose one filling method based on user input*/
void choose_fill_strategy(Model_F &fasta_in,ModelUnit &pdb_in,map<string, map<string,vector<double> > >&ideal_rna,vector<pair<double,vector<size_t> > >&bp_vec,vector<vector<size_t> > &res_str_vec,const int option=6){
    while(can_handle(fasta_in,pdb_in,option)){
        if(option==0) fill_gap_whole(fasta_in,pdb_in,option);
        else if(option==1) extrapolate(fasta_in,pdb_in);
        else if(option==2) fill_gap_new(fasta_in,pdb_in,option);
        else if(option==3) fill_gap(fasta_in,pdb_in);
        else if(option==4) fill_gap_new(fasta_in,pdb_in,option);
        else if(option==5) fill_gap_new(fasta_in,pdb_in,option);
        else if(option==6){
            fill_arc(fasta_in,pdb_in);
            MissingRNAatom(pdb_in,ideal_rna,bp_vec,7);
            standardize_residue_order(pdb_in);
            adjust_position(pdb_in,res_str_vec,bp_vec,ideal_rna);
            MissingRNAatom(pdb_in,ideal_rna,bp_vec,6);
            adjust_position(pdb_in,res_str_vec,bp_vec,ideal_rna);
            MissingRNAatom(pdb_in,ideal_rna,bp_vec,5);
        }
        else if(option==7){
            fill_SRS(fasta_in,pdb_in);
            MissingRNAatom(pdb_in,ideal_rna,bp_vec,7);
            standardize_residue_order(pdb_in);
            adjust_position(pdb_in,res_str_vec,bp_vec,ideal_rna);
            MissingRNAatom(pdb_in,ideal_rna,bp_vec,6);
            adjust_position(pdb_in,res_str_vec,bp_vec,ideal_rna);
            MissingRNAatom(pdb_in,ideal_rna,bp_vec,5);
        }
    }
    
}

#endif
