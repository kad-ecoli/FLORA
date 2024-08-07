/*To */

#ifndef SSS_HPP
#define SSS_HPP
#include <iostream>
#include <map>
#include <vector>
#include <stack>
#include <fstream>
#include "PDBParser.hpp"
#include "PDBFiller.hpp"
#include "Superpose.hpp"
#include "MissingRNAatom.hpp"
#include "IdealRNA.hpp"

using namespace std;

bool quick_superpose(ResidueUnit &residue,map<string, map<string,vector<double> > > &ideal_rna){
    vector<double> tmp(3,0);
    vector<vector<double> > xyz_list1;
    vector<vector<double> > xyz_list2;
    vector<vector<double> > RotMatix;
    vector<double> TranVect;
    //cout<<residue.resn<<endl;
    string key=residue.resn;
    size_t atomNum=residue.atoms.size();
    xyz_list1.assign(atomNum,tmp);
    xyz_list2.assign(atomNum,tmp);
    for(size_t i=0;i<residue.atoms.size();i++){
        string atom_name=residue.atoms[i].name;
        xyz_list1[i][0]=ideal_rna[key][atom_name][0];
        xyz_list1[i][1]=ideal_rna[key][atom_name][1];
        xyz_list1[i][2]=ideal_rna[key][atom_name][2];

        xyz_list2[i][0]=residue.atoms[i].xyz[0];
        xyz_list2[i][1]=residue.atoms[i].xyz[1];
        xyz_list2[i][2]=residue.atoms[i].xyz[2];
    }
    atomNum=countUniqAtom(xyz_list2);
    //cout<<atomNum<<endl;
    if (atomNum<xyz_list2.size()) cerr<<"duplicated coordinate when reconstructing"<<endl;
    if (atomNum<3)
    {
        for (int j=0;j<xyz_list1.size();j++)
        {
            xyz_list1[j].clear();
            xyz_list1[j].shrink_to_fit();
            xyz_list2[j].clear();
            xyz_list2[j].shrink_to_fit();
        }
        xyz_list1.clear();
        xyz_list1.shrink_to_fit();
        xyz_list2.clear();
        xyz_list2.shrink_to_fit();
        key.clear();
        return false;
    }
    RotateCoor(xyz_list1,xyz_list2, RotMatix, TranVect);

    AtomUnit atom;
    atom.movable=1;
    atom.name="";
    atom.xyz.assign(3,0);
    vector<AtomUnit>().swap(residue.atoms);
    residue.atoms.assign(ideal_rna[key].size(),atom);
    int a=0;
    for(auto iter:ideal_rna[key]){
        residue.atoms[a].name=iter.first;
        ChangeCoor(iter.second,RotMatix,TranVect,residue.atoms[a].xyz);
        a++;
    }
    //cout<<"good"<<endl;
    standardize_atom_order(residue);
    //cout<<"goodagain"<<endl;

    for (int j=0;j<xyz_list1.size();j++)
    {
        xyz_list1[j].clear();
        xyz_list1[j].shrink_to_fit();
        xyz_list2[j].clear();
        xyz_list2[j].shrink_to_fit();
    }
    xyz_list1.clear();
    xyz_list1.shrink_to_fit();
    xyz_list2.clear();
    xyz_list2.shrink_to_fit();
    TranVect.clear();
    TranVect.shrink_to_fit();
    for (int j=0;j<3;j++) RotMatix[j].clear();
    for (int j=0;j<3;j++) RotMatix[j].shrink_to_fit();
    RotMatix.clear();
    RotMatix.shrink_to_fit();
    key.clear();
    key.shrink_to_fit();
    //cout<<"still_good"<<endl;
    return true;
}

Model_F read_dbn(char* filename,const int option=0){
    Model_F dbn;
    ifstream fp;
    fp.open(filename,ios::in);
    if(option==0){
    string line;
    int chainI=-1;
    while(fp.good()){
        getline(fp,line);
        if(!fp.good()) break;
        if(line[0]=='>'){
            chainI++;
            continue;
        }
        if(dbn.chains_F.size()==chainI) dbn.chains_F.push_back(line);
        else if(dbn.chains_F.size()>chainI) dbn.chains_F[chainI]+=line;
    }
    }
    else if(option==1){
        string line;
        getline(fp,line);
        dbn.chains_F.push_back(line);
    }
    fp.close();
    //cout<<dbn.chains_F[0]<<endl;
    return dbn;
}

vector<int> pair_info(string dot_sequence){
    map<int,int> pairs;
    stack<int> agbk;
    stack<int> robk;
    stack<int> sqbk;
    stack<int> cubk;
    //cout<<"Hello"<<endl;
    string temp=dot_sequence;
    dot_sequence="";
    for(int i=0;i<temp.size();i++){
        if(temp[i]!='&') dot_sequence+=temp[i];
    }
    for(int i=0;i<dot_sequence.size();i++){
        if(dot_sequence[i]=='<') agbk.push(i+1);
        else if(dot_sequence[i]=='(') robk.push(i+1);
        else if(dot_sequence[i]=='[') sqbk.push(i+1);
        else if(dot_sequence[i]=='{') cubk.push(i+1);
        else if(dot_sequence[i]=='>'){
            pairs.insert(pair<int,int>(i+1,agbk.top()));
            pairs.insert(pair<int,int>(agbk.top(),i+1));
            agbk.pop();
        }
        else if(dot_sequence[i]==')'){
            pairs.insert(pair<int,int>(i+1,robk.top()));
            pairs.insert(pair<int,int>(robk.top(),i+1));
            robk.pop();
        }
        else if(dot_sequence[i]==']'){
            pairs.insert(pair<int,int>(i+1,sqbk.top()));
            pairs.insert(pair<int,int>(sqbk.top(),i+1));
            sqbk.pop();
        }
        else if(dot_sequence[i]=='}'){
            pairs.insert(pair<int,int>(i+1,cubk.top()));
            pairs.insert(pair<int,int>(cubk.top(),i+1));
            cubk.pop();
        }
    }
    vector<int> basepair;
    for(int i=1;i<=dot_sequence.size();i++){
        if(pairs.find(i)!=pairs.end()) basepair.push_back(pairs[i]);
        else basepair.push_back(0);
    }
    map<int,int>().swap(pairs);
    //cout<<"pair_info is good!"<<endl;
    return basepair;
}

bool fill_single_pair(int solid_resi_num, int empty_resi_num, string empty_resi_name, 
                      ChainUnit &chain_in,map<string, map<string,vector<double> > > &ideal_rna){
    vector<double> tmp(3,0);
    vector<vector<double> > xyz_list1;
    vector<vector<double> > xyz_list2;
    vector<vector<double> > RotMatix;
    vector<double> TranVect;
    string key1="";
    int num=0;
    for(num=0;num<chain_in.residues.size();num++){
        if(chain_in.residues[num].resi==solid_resi_num){
            key1=chain_in.residues[num].resn;
            break;
        } 

    }
    if(key1==""){
        cout<<solid_resi_num<<endl;
        return false;
    } 
    //cout<<solid_resi_num<<" "<<num<<" "<<chain_in.residues[num].resi<<endl;
    string key2=empty_resi_name;

    if(key1=="  G"&&key2=="  U"){
        key1="U G";key2="G U";
    }

    if(key1=="  U"&&key2=="  G"){
        key1="G U";key2="U G";
    }

    int a,a1,a2;
    int atomNum=0;
    atomNum=chain_in.residues[num].atoms.size();
    xyz_list1.assign(atomNum,tmp);
    xyz_list2.assign(atomNum,tmp);
    vector<double>().swap(tmp);
    a=0;
    for (a1=0;a1<chain_in.residues[num].atoms.size();a1++)
    {
        xyz_list1[a][0]=ideal_rna[key1][chain_in.residues[num].atoms[a1].name][0];
        xyz_list1[a][1]=ideal_rna[key1][chain_in.residues[num].atoms[a1].name][1];
        xyz_list1[a][2]=ideal_rna[key1][chain_in.residues[num].atoms[a1].name][2];

        xyz_list2[a][0]=chain_in.residues[num].atoms[a1].xyz[0];
        xyz_list2[a][1]=chain_in.residues[num].atoms[a1].xyz[1];
        xyz_list2[a][2]=chain_in.residues[num].atoms[a1].xyz[2];
        a++;
    }
    atomNum=countUniqAtom(xyz_list2);
    if (atomNum<xyz_list2.size()) cerr<<"duplicated coordinate when reconstructing"<<endl;
    if (atomNum<3)
    {
        for (a=0;a<xyz_list1.size();a++)
        {
            xyz_list1.clear();
            xyz_list1.shrink_to_fit();
            xyz_list2.clear();
            xyz_list2.shrink_to_fit();
        }
        xyz_list1.clear();
        xyz_list1.shrink_to_fit();
        xyz_list2.clear();
        xyz_list2.shrink_to_fit();
        key1.clear();
        key2.clear();
        return false;
    }
    RotateCoor(xyz_list1,xyz_list2, RotMatix, TranVect);
    ResidueUnit New_Resi;
    AtomUnit fill_atom;
    New_Resi.het=false;
    New_Resi.icode='1';
    New_Resi.resi=empty_resi_num;
    New_Resi.resn=empty_resi_name;
    fill_atom.name="";
    fill_atom.movable=1;
    fill_atom.xyz.assign(3,0);
    //cout<<ideal_rna[key2].size()<<endl;
    New_Resi.atoms.assign(ideal_rna[key2].size(),fill_atom);
    a=0;
    for (auto iter:ideal_rna[key2]){
        //cout<<iter.first<<'\t'<<iter.second[0]<<'\t'<<iter.second[1]<<'\t'<<iter.second[2]<<endl;
        New_Resi.atoms[a].name=iter.first;
        ChangeCoor(iter.second,RotMatix,TranVect,New_Resi.atoms[a].xyz);
        a++;
    }
    standardize_atom_order(New_Resi);
    
    /* clean up */
    for (a=0;a<xyz_list1.size();a++)
    {
        xyz_list1[a].clear();
        xyz_list1[a].shrink_to_fit();
        xyz_list2[a].clear();
        xyz_list2[a].shrink_to_fit();
    }
    xyz_list1.clear();
    xyz_list1.shrink_to_fit();
    xyz_list2.clear();
    xyz_list2.shrink_to_fit();
    TranVect.clear();
    TranVect.shrink_to_fit();
    for (a=0;a<3;a++) RotMatix[a].clear();
    for (a=0;a<3;a++) RotMatix[a].shrink_to_fit();
    RotMatix.clear();
    RotMatix.shrink_to_fit();
    key1.clear();
    key2.clear();
    chain_in.residues.push_back(New_Resi);
    standardize_residue_order(chain_in);
    vector<AtomUnit>().swap(New_Resi.atoms);
    //cout<<"put_single_pair is good!"<<endl;
    return true;
}

void single_pairs(vector<int> &pairif, sizet_sizet_vc &gapif, string &sequence, ChainUnit &chain_in, map<string, map<string,vector<double> > > &ideal_rna){
    map<int,pair<int,string> > pairs;
    map<int,pair<int,string> > pairs_s;
    for(int j=0;j<gapif.size();j++){
        int key;
        string value="   ";
        for(int k=gapif[j].first;k<gapif[j].second;k++){
            key=pairif[k-1];
            if(key){
                value[2]=toupper(sequence[k-1]);
                pairs.insert(pair<int,pair<int,string> >(key,pair<int,string>(k,value)));
            }
        }
    }
    for(auto iter:pairs){
        int in=0;
        for(int j=0;j<gapif.size();j++){
            if(iter.first>=gapif[j].first&&iter.first<gapif[j].second) in++;
        }
        if(in==0) pairs_s.insert(iter);
    }
    //cout<<pairs_s.size()<<endl;
    map<int,pair<int,string> >().swap(pairs);
    //for(auto iter:pairs_s) cout<<iter.first<<"\t"<<iter.second.first<<"\t"<<iter.second.second<<endl;
    for(auto iter:pairs_s) fill_single_pair(iter.first,iter.second.first,iter.second.second,chain_in,ideal_rna);
    map<int,pair<int,string> >().swap(pairs_s);
    //cout<<"single_pair is good!"<<endl;
}

void extrapolate_single_pairs(string &sequence, ChainUnit &chain_in,int start_resi,int helix_length, vector<bool> &existence){//if it is extended from both trans side, there is no need of input end_resi
    cout<<"start to extrapolate:"<<start_resi<<"\t"<<helix_length<<endl;
    vector<int> helix_gap;//start+length+f_num+r_num+possible_to_extrapolate
    vector<vector<int> > helix_gaps;
    helix_gap.assign(5,0);
    bool sw=true;
    int f=0;
    int l=0;
    for(int i=0;i<existence.size();i++){
        if(existence[i]&&!l){
            f++;
        }else if(!existence[i]&&!l){
            if(helix_gap[1]){
                helix_gap[3]=f;
                helix_gaps.push_back(helix_gap);
            }
            helix_gap[0]=i;
            helix_gap[2]=f;
            f=0;
            l++;
        }else if(!existence[i]&&l){
            l++;
        }else if(existence[i]&&l){
            helix_gap[1]=l;
            l=0;
            f++;
        }
    }
    if(l) helix_gap[1]=l;
    helix_gap[3]=f;
    helix_gaps.push_back(helix_gap);
    /*if(l){
        helix_gap[1]=l;
        helix_gap[3]=f;
        helix_gaps.push_back(helix_gap);
    }
    if(f){
        helix_gap[3]=f;
        helix_gaps.push_back(helix_gap);
    }*/
    cout<<"helix_gaps_size:"<<helix_gaps.size()<<endl;
    for(int i=0;i<helix_gaps.size();i++){
        cout<<"helix_gaps:"<<helix_gaps[i][0]+start_resi<<"\t"<<helix_gaps[i][1]<<"\t"<<helix_gaps[i][2]<<"\t"<<helix_gaps[i][3]<<"\t"<<helix_gaps[i][4]<<endl;
        if(helix_gaps[i][1]>0&&(helix_gaps[i][2]>1||helix_gaps[i][3]>1)){
            if(helix_gaps[i][2]>1){
                vector<ResidueUnit> resi_fs;
                int num_f=4;
                if(num_f>helix_gaps[i][2]) num_f=helix_gaps[i][2];
                for(int r=0;r<chain_in.residues.size();r++){
                    for(int fs=0;fs<num_f;fs++){
                        if(chain_in.residues[r].resi==start_resi+helix_gaps[i][0]-num_f+fs) resi_fs.push_back(chain_in.residues[r]);
                        if(chain_in.residues[r].resi==start_resi+helix_gaps[i][0]-num_f+fs) cout<<chain_in.residues[r].resi<<endl;
                    }
                }
                cout<<"reference_resi_num:"<<resi_fs.size()<<endl;
                vector<ResidueUnit> new_resi;
                ResidueUnit resi_temp;
                resi_temp.het=false;
                resi_temp.icode='1';
                resi_temp.resn="   ";
                AtomUnit atom_temp;
                atom_temp.movable=1;
                atom_temp.xyz.assign(3,0);
                resi_temp.atoms.assign(12,atom_temp);
                vector<vector<double> > xyz_listt;
                vector<vector<double> > xyz_listr;
                xyz_listt.assign(12*(num_f-1),atom_temp.xyz);
                xyz_listr.assign(12*(num_f-1),atom_temp.xyz);
                vector<vector<double> > RotMatix;
                vector<double> TranVect;
                for(int g=0;g<helix_gaps[i][1];g++){
                    resi_temp.resi=resi_fs[resi_fs.size()-1].resi+1;
                    resi_temp.resn[2]=toupper(sequence[resi_temp.resi-1]);
                    for(int j=0;j<num_f-1;j++){
                        for(int k=0;k<12;k++){
                            for(int l=0;l<3;l++){
                                xyz_listt[12*j+k][l]=resi_fs[resi_fs.size()-num_f+j].atoms[k].xyz[l];
                                xyz_listr[12*j+k][l]=resi_fs[resi_fs.size()-num_f+1+j].atoms[k].xyz[l];
                            }
                        }
                    }
                    
                    RotateCoor(xyz_listt,xyz_listr, RotMatix, TranVect);
                    for(int k=0;k<12;k++){
                        resi_temp.atoms[k].name=resi_fs[resi_fs.size()-1].atoms[k].name;
                        ChangeCoor(resi_fs[resi_fs.size()-1].atoms[k].xyz,RotMatix,TranVect,resi_temp.atoms[k].xyz);
                    }
                    resi_fs.push_back(resi_temp);
                    new_resi.push_back(resi_temp);
                }
                cout<<"gap_length and new_resi size:"<<helix_gaps[i][1]<<"\t"<<new_resi.size()<<endl;
                fill_new_residue(chain_in,new_resi);
                vector<vector<double> >().swap(xyz_listt);
                vector<vector<double> >().swap(xyz_listr);
                vector<ResidueUnit>().swap(resi_fs);
                vector<AtomUnit>().swap(resi_temp.atoms);
                vector<ResidueUnit>().swap(new_resi);
            }else if(helix_gaps[i][3]>1){
                vector<ResidueUnit> resi_rs;
                int num_r=4;
                if(num_r>helix_gaps[i][3]) num_r=helix_gaps[i][3];
                for(int r=chain_in.residues.size()-1;r>=0;r--){
                    for(int rs=0;rs<num_r;rs++){
                        if(chain_in.residues[r].resi==start_resi+helix_gaps[i][0]+helix_gaps[i][1]+num_r-1-rs) resi_rs.push_back(chain_in.residues[r]);
                        if(chain_in.residues[r].resi==start_resi+helix_gaps[i][0]+helix_gaps[i][1]+num_r-1-rs) cout<<chain_in.residues[r].resi<<endl;
                    }
                }//the order of residue is decided by r++/r--
                cout<<"reference_resi_num:"<<resi_rs.size()<<endl;
                vector<ResidueUnit> new_resi;
                ResidueUnit resi_temp;
                resi_temp.het=false;
                resi_temp.icode='1';
                resi_temp.resn="   ";
                AtomUnit atom_temp;
                atom_temp.movable=1;
                atom_temp.xyz.assign(3,0);
                resi_temp.atoms.assign(12,atom_temp);
                vector<vector<double> > xyz_listt;
                vector<vector<double> > xyz_listr;
                xyz_listt.assign(12*(num_r-1),atom_temp.xyz);
                xyz_listr.assign(12*(num_r-1),atom_temp.xyz);
                vector<vector<double> > RotMatix;
                vector<double> TranVect;
                for(int g=0;g<helix_gaps[i][1];g++){
                    resi_temp.resi=resi_rs[resi_rs.size()-1].resi-1;
                    cout<<"resi_temp.resi:"<<resi_temp.resi<<endl;
                    resi_temp.resn[2]=toupper(sequence[resi_temp.resi-1]);
                    for(int j=0;j<num_r-1;j++){
                        for(int k=0;k<12;k++){
                            for(int l=0;l<3;l++){
                                xyz_listt[12*j+k][l]=resi_rs[resi_rs.size()-num_r+j].atoms[k].xyz[l];
                                xyz_listr[12*j+k][l]=resi_rs[resi_rs.size()-num_r+1+j].atoms[k].xyz[l];
                            }
                        }
                    }
                    
                    RotateCoor(xyz_listt,xyz_listr, RotMatix, TranVect);
                    for(int k=0;k<12;k++){
                        resi_temp.atoms[k].name=resi_rs[resi_rs.size()-1].atoms[k].name;
                        ChangeCoor(resi_rs[resi_rs.size()-1].atoms[k].xyz,RotMatix,TranVect,resi_temp.atoms[k].xyz);
                    }
                    resi_rs.push_back(resi_temp);
                    new_resi.push_back(resi_temp);
                }
                cout<<"gap_length and new_resi size:"<<helix_gaps[i][1]<<"\t"<<new_resi.size()<<endl;
                fill_new_residue(chain_in,new_resi);
                vector<vector<double> >().swap(xyz_listt);
                vector<vector<double> >().swap(xyz_listr);
                vector<ResidueUnit>().swap(resi_rs);
                vector<AtomUnit>().swap(resi_temp.atoms);
                vector<ResidueUnit>().swap(new_resi);
            }
            standardize_residue_order(chain_in);
        }
    }
}

void extrapolate_single_pairs(vector<int> &pairif, string &sequence, ChainUnit &chain_in){
    map<int,pair<int,int> > helix_info;
    vector<int> helix_info_temp;
    helix_info_temp.assign(3,0);
    for(int i=0;i<pairif.size();i++){
        if(pairif[i]==0) continue;
        if(i+1>helix_info_temp[0]+helix_info_temp[2]){
            helix_info[helix_info_temp[0] ]=pair<int,int>(helix_info_temp[1],helix_info_temp[2]);
            helix_info_temp[0]=i+1;
            helix_info_temp[1]=pairif[i];
            helix_info_temp[2]=1;
        }else
        if(i+1==helix_info_temp[0]+helix_info_temp[2]){
            helix_info_temp[1]--;
            helix_info_temp[2]++;
        }
        cout<<helix_info_temp[0]<<"\t"<<helix_info_temp[1]<<"\t"<<helix_info_temp[2]<<endl;
    }
    helix_info[helix_info_temp[0] ]=pair<int,int>(helix_info_temp[1],helix_info_temp[2]);
    map<int,pair<int,int> > helix_info_s;
    for(auto iter:helix_info){
        if(!iter.first) continue;
        if(helix_info_s.find(iter.second.first-iter.second.second+1)==helix_info_s.end()) helix_info_s.insert(iter);
    }
    cout<<helix_info.size()-1<<"\t"<<helix_info_s.size()<<endl;
    map<int,pair<int,int> >().swap(helix_info_s);
    vector<vector<bool> > helix_info_d;
    vector<bool> helix_info_s_d;
    for(auto iter:helix_info){//_s
        vector<bool>().swap(helix_info_s_d);
        helix_info_s_d.assign(iter.second.second,false);
        for(int j=iter.first;j<iter.first+iter.second.second;j++){
            for(int k=0;k<chain_in.residues.size();k++){
                if(chain_in.residues[k].resi==j){
                    helix_info_s_d[j-iter.first]=true;
                    break;
                }
            }
        }
        helix_info_d.push_back(helix_info_s_d);
    }
    vector<bool>().swap(helix_info_s_d);
    cout<<helix_info_d.size()<<"\t"<<helix_info_s.size()<<endl;
    int i=0;
    for(auto iter:helix_info){
        extrapolate_single_pairs(sequence,chain_in,iter.first,iter.second.second,helix_info_d[i]);
        i++;
    }
    vector<vector<bool> >().swap(helix_info_d);
    map<int,pair<int,int> >().swap(helix_info);
}

void fill_double_pair_backbone(map<int,int> &pairs_s,ChainUnit &chain_in,map<string, map<string,vector<double> > > &ideal_rna){
    vector<vector<int> > helix_info;
    vector<int> helix_info_temp;
    helix_info_temp.assign(3,0);
    for(auto iter:pairs_s){
        if(iter.first>helix_info_temp[0]+helix_info_temp[2]){
            helix_info.push_back(helix_info_temp);
            helix_info_temp[0]=iter.first;
            helix_info_temp[1]=iter.second;
            helix_info_temp[2]=1;
        }else
        if(iter.first==helix_info_temp[0]+helix_info_temp[2]){
            helix_info_temp[1]--;
            helix_info_temp[2]++;
        }
        cout<<helix_info_temp[0]<<"\t"<<helix_info_temp[1]<<"\t"<<helix_info_temp[2]<<endl;
    }
    helix_info.push_back(helix_info_temp);
    int max_helix_length=48;
    for(int i=1;i<helix_info.size();i++){
        if(helix_info[i][2]>max_helix_length){
            helix_info_temp[0]=helix_info[i][0]+max_helix_length;
            helix_info_temp[1]=helix_info[i][1];
            helix_info_temp[2]=helix_info[i][2]-max_helix_length;
            helix_info[i][1]+=helix_info[i][2];
            helix_info[i][1]-=max_helix_length;
            helix_info[i][2]=max_helix_length;
            helix_info.push_back(helix_info_temp);
        }
    }
    //cout<<helix_info.size()<<endl;
    vector<ResidueUnit> new_resi;
    ResidueUnit resi_temp;
    AtomUnit atom_temp;
    atom_temp.xyz.assign(3,0);
    for(int i=1;i<helix_info.size();i++){
        vector<vector<double> > xyz_list1;
        //cout<<helix_info[i][2]<<"\t"<<xyz_list1.size()<<endl;
        vector<vector<double> > xyz_list2;
        vector<vector<double> > RotMatix;
        vector<double> TranVect;
        for(int r=0;r<chain_in.residues.size();r++){
            if(chain_in.residues[r].resi>=helix_info[i][0]&&chain_in.residues[r].resi<helix_info[i][0]+helix_info[i][2]){
                for(int a=0;a<12;a++){
                    xyz_list2.push_back(chain_in.residues[r].atoms[a].xyz);
                }
            }
            if(chain_in.residues[r].resi>=helix_info[i][1]&&chain_in.residues[r].resi<helix_info[i][1]+helix_info[i][2]){
                for(int a=0;a<12;a++){
                    xyz_list2.push_back(chain_in.residues[r].atoms[a].xyz);
                }
            }
        }
        if(xyz_list2.size()<12*2*helix_info[i][2]) continue;
        xyz_list1=ideal_RNA_backbone(helix_info[i][2]);
        //cout<<xyz_list2.size()<<endl;
        RotateCoor(xyz_list1,xyz_list2, RotMatix, TranVect);
        int atom_num=0;
        for(int r=0;r<chain_in.residues.size();r++){
            if(chain_in.residues[r].resi>=helix_info[i][0]&&chain_in.residues[r].resi<helix_info[i][0]+helix_info[i][2]){
                resi_temp.het=chain_in.residues[r].het;
                resi_temp.icode='1';
                resi_temp.resi=chain_in.residues[r].resi;
                resi_temp.resn=chain_in.residues[r].resn;
                atom_temp.movable=1;
                resi_temp.atoms.assign(12,atom_temp);
                for(int a=0;a<12;a++){
                    resi_temp.atoms[a].name=chain_in.residues[r].atoms[a].name;
                    ChangeCoor(xyz_list1[atom_num],RotMatix,TranVect,resi_temp.atoms[a].xyz);
                    atom_num++;
                }
                if(chain_in.residues[r].icode!='1') new_resi.push_back(resi_temp);
                vector<AtomUnit>().swap(resi_temp.atoms);
            }
            if(chain_in.residues[r].resi>=helix_info[i][1]&&chain_in.residues[r].resi<helix_info[i][1]+helix_info[i][2]){
                resi_temp.het=chain_in.residues[r].het;
                resi_temp.icode='1';
                resi_temp.resi=chain_in.residues[r].resi;
                resi_temp.resn=chain_in.residues[r].resn;
                atom_temp.movable=1;
                resi_temp.atoms.assign(12,atom_temp);
                for(int a=0;a<12;a++){
                    resi_temp.atoms[a].name=chain_in.residues[r].atoms[a].name;
                    ChangeCoor(xyz_list1[atom_num],RotMatix,TranVect,resi_temp.atoms[a].xyz);
                    atom_num++;
                }
                if(chain_in.residues[r].icode!='1') new_resi.push_back(resi_temp);
                vector<AtomUnit>().swap(resi_temp.atoms);
            }
        }
        xyz_list1.clear();
        xyz_list1.shrink_to_fit();
        xyz_list2.clear();
        xyz_list2.shrink_to_fit();
        TranVect.clear();
        TranVect.shrink_to_fit();
        for (int a=0;a<3;a++) RotMatix[a].clear();
        for (int a=0;a<3;a++) RotMatix[a].shrink_to_fit();
        RotMatix.clear();
        RotMatix.shrink_to_fit();
    }
    delete_moveable_residue(chain_in);
    fill_new_residue(chain_in,new_resi);
    vector<ResidueUnit>().swap(new_resi);
    standardize_residue_order(chain_in);
}
bool fill_double_pair(int first_resi_num,string first_resi_name, int second_resi_num, string second_resi_name, 
    ChainUnit &chain_in, string sequence, map<string, map<string,vector<double> > > &ideal_rna){
    standardize_residue_order(chain_in);
    vector<double> tmp(3,0);
    vector<vector<double> > xyz_list1;
    vector<vector<double> > xyz_list2;
    vector<vector<double> > RotMatix;
    vector<double> TranVect;
    string key1=first_resi_name;
    int num1=0;
    for(num1=0;num1<chain_in.residues.size();num1++){
        if(chain_in.residues[num1].resi==first_resi_num) break;
    }
    //cout<<first_resi_num<<" "<<num1<<" "<<chain_in.residues[num1].resi<<endl;
    
    string key2=second_resi_name;
    int num2=0;
    for(num2=0;num2<chain_in.residues.size();num2++){
        if(chain_in.residues[num2].resi==second_resi_num) break;
    }
    //cout<<second_resi_num<<" "<<num2<<" "<<chain_in.residues[num2].resi<<endl;

    size_t a,a1,a2;
    size_t atomNum=0;
    atomNum=chain_in.residues[num1].atoms.size()+chain_in.residues[num2].atoms.size();
    xyz_list1.assign(atomNum,tmp);
    xyz_list2.assign(atomNum,tmp);
    a=0;
    for (a1=0;a1<chain_in.residues[num1].atoms.size();a1++)
    {
        xyz_list1[a][0]=ideal_rna[key1][chain_in.residues[num1].atoms[a1].name][0];
        xyz_list1[a][1]=ideal_rna[key1][chain_in.residues[num1].atoms[a1].name][1];
        xyz_list1[a][2]=ideal_rna[key1][chain_in.residues[num1].atoms[a1].name][2];

        xyz_list2[a][0]=chain_in.residues[num1].atoms[a1].xyz[0];
        xyz_list2[a][1]=chain_in.residues[num1].atoms[a1].xyz[1];
        xyz_list2[a][2]=chain_in.residues[num1].atoms[a1].xyz[2];
        //cout<<chain_in.residues[num1].atoms[a1].name<<'\t'
        //    <<xyz_list1[a][0]<<'\t'<<xyz_list1[a][1]<<'\t'<<xyz_list1[a][2]<<'\t'
        //    <<xyz_list2[a][0]<<'\t'<<xyz_list2[a][1]<<'\t'<<xyz_list2[a][2]<<endl;
        a++;
    }
    //cout<<num1<<'\t'<<num2<<endl;
    for (a2=0;a2<chain_in.residues[num2].atoms.size();a2++)
    {
        xyz_list1[a][0]=ideal_rna[key2][chain_in.residues[num2].atoms[a2].name][0];
        xyz_list1[a][1]=ideal_rna[key2][chain_in.residues[num2].atoms[a2].name][1];
        xyz_list1[a][2]=ideal_rna[key2][chain_in.residues[num2].atoms[a2].name][2];

        xyz_list2[a][0]=chain_in.residues[num2].atoms[a2].xyz[0];
        xyz_list2[a][1]=chain_in.residues[num2].atoms[a2].xyz[1];
        xyz_list2[a][2]=chain_in.residues[num2].atoms[a2].xyz[2];
        //cout<<chain_in.residues[num2].atoms[a2].name<<'\t'
        //    <<xyz_list1[a][0]<<'\t'<<xyz_list1[a][1]<<'\t'<<xyz_list1[a][2]<<'\t'
        //    <<xyz_list2[a][0]<<'\t'<<xyz_list2[a][1]<<'\t'<<xyz_list2[a][2]<<endl;
        a++;
    }
    atomNum=countUniqAtom(xyz_list2);
    if (atomNum<xyz_list2.size()) cerr<<"duplicated coordinate when reconstructing"<<endl;
    if (atomNum<3)
    {
        for (a=0;a<xyz_list1.size();a++)
        {
            xyz_list1[a].clear();
            xyz_list1[a].shrink_to_fit();
            xyz_list2[a].clear();
            xyz_list2[a].shrink_to_fit();
        }
        xyz_list1.clear();
        xyz_list1.shrink_to_fit();
        xyz_list2.clear();
        xyz_list2.shrink_to_fit();
        key1.clear();
        key1.shrink_to_fit();
        key2.clear();
        key2.shrink_to_fit();
        return false;
    }
    RotateCoor(xyz_list1,xyz_list2, RotMatix, TranVect);
    //cout<<RotMatix[0][0]<<'\t'<<RotMatix[0][1]<<'\t'<<RotMatix[0][2]<<'\t'<<TranVect[0]<<endl
    //    <<RotMatix[1][0]<<'\t'<<RotMatix[1][1]<<'\t'<<RotMatix[1][2]<<'\t'<<TranVect[1]<<endl
    //    <<RotMatix[2][0]<<'\t'<<RotMatix[2][1]<<'\t'<<RotMatix[2][2]<<'\t'<<TranVect[2]<<endl;
    
    for (int i=0;i<chain_in.residues[num1].atoms.size();i++){
        ChangeCoor(ideal_rna[key1][chain_in.residues[num1].atoms[i].name],RotMatix,TranVect,chain_in.residues[num1].atoms[i].xyz);
    }

    for(int i=0;i<chain_in.residues[num2].atoms.size();i++){
        ChangeCoor(ideal_rna[key2][chain_in.residues[num2].atoms[i].name],RotMatix,TranVect,chain_in.residues[num2].atoms[i].xyz);
    }
    
    
    for (a=0;a<xyz_list1.size();a++)
    {
        xyz_list1[a].clear();
        xyz_list1[a].shrink_to_fit();
        xyz_list2[a].clear();
        xyz_list2[a].shrink_to_fit();
    }
    xyz_list1.clear();
    xyz_list1.shrink_to_fit();
    xyz_list2.clear();
    xyz_list2.shrink_to_fit();
    TranVect.clear();
    TranVect.shrink_to_fit();
    for (a=0;a<3;a++) RotMatix[a].clear();
    for (a=0;a<3;a++) RotMatix[a].shrink_to_fit();
    RotMatix.clear();
    RotMatix.shrink_to_fit();
    key1.clear();
    key1.shrink_to_fit();
    key2.clear();
    key2.shrink_to_fit();
    standardize_residue_order(chain_in);

    //cout<<num1<<" "<<num2<<endl;

    int m=0;
    int n=0;
    if(num1<num2){
        m=num1;
        n=num2;
    }else{
        m=num2;
        n=num1;
    }
    vector<ResidueUnit> new_resi;
    if(m>0&&chain_in.residues[m].resi==chain_in.residues[m-1].resi+1&&
       n<chain_in.residues.size()-1&&chain_in.residues[n+1].resi==chain_in.residues[n].resi+1
       &&Points2Distance(chain_in.residues[m].atoms[0].xyz,chain_in.residues[n].atoms[0].xyz)<5.4*(chain_in.residues[n].resi-chain_in.residues[m].resi)){
        vector<vector<double> > adjacent_atoms;
        adjacent_atoms.push_back(chain_in.residues[m-1].atoms[5].xyz);
        adjacent_atoms.push_back(chain_in.residues[m].atoms[5].xyz);
            adjacent_atoms.push_back(chain_in.residues[n].atoms[5].xyz);
            adjacent_atoms.push_back(chain_in.residues[n+1].atoms[5].xyz);
            //cout<<chain_in.residues[m].resi<<endl;
            RNAArc(Points2Distance(adjacent_atoms[1],adjacent_atoms[2]),chain_in.residues[n].resi-chain_in.residues[m].resi,chain_in.residues[m].resi,sequence,adjacent_atoms,new_resi);
            vector<vector<double> >().swap(adjacent_atoms);
       }
    else fill_linear(sequence,chain_in.residues[m],chain_in.residues[n],new_resi);
    for(int i=0;i<new_resi.size();i++){
        quick_superpose(new_resi[i],ideal_rna);
        standardize_atom_order(new_resi[i]);
    }
    replace_new_residue(chain_in,new_resi);
    vector<ResidueUnit>().swap(new_resi);
    
    //cout<<"put_double_pair is good!"<<endl;
    return true;
}

void double_pairs(vector<int> &pairif, sizet_sizet_vc &gapif, string &sequence, ChainUnit &chain_in, map<string, map<string,vector<double> > > &ideal_rna){
    map<int,int> pairs;
    map<int,int> pairs_s;
    //cout<<pairif.size()<<endl;
    //cout<<gapif.size()<<endl;
    for(int i=0;i<pairif.size();i++){
        //cout<<i+1<<"\t"<<pairif[i]<<endl;
    }
    for(int j=0;j<gapif.size();j++){
        int key;
        string value="   ";
        for(int k=gapif[j].first;k<gapif[j].second;k++){
            key=pairif[k-1];
            //cout<<k<<" "<<key<<endl;
            if(key&&pairs.find(k)==pairs.end()){
                pairs.insert(pair<int,int>(key,k));
            }
        }
    }
    //cout<<pairs.size()<<endl;
    for(auto iter:pairs){
        int in=0;
        for(int j=0;j<gapif.size();j++){
            if(iter.first>=gapif[j].first&&iter.first<gapif[j].second) in++;
            if(iter.second>=gapif[j].first&&iter.second<gapif[j].second) in++;
        }
        //cout<<in<<endl;
        if(in==2) pairs_s.insert(pair<int,int>(iter.second,iter.first));
    }
    fill_double_pair_backbone(pairs_s,chain_in,ideal_rna);
    /*string value1="   ";
    string value2="   ";
    for(auto iter:pairs_s){
        cout<<iter.first<<"\t"<<iter.second<<endl;
        value1[2]=toupper(sequence[iter.first-1]);
        value2[2]=toupper(sequence[iter.second-1]);
        fill_double_pair(iter.first,value1,iter.second,value2,chain_in,sequence,ideal_rna);
        string name=to_string(iter.first)+"_"+to_string(iter.second)+".pdb";
        write_pdb_structure(name.c_str(),chain_in);
    }*/
    map<int,int>().swap(pairs_s);
    //cout<<"double_pair is good!"<<endl;
}

vector<sizet_sizet_vc> put_single_pairs(char* filename, Model_F &fasta_in, ModelUnit &pdb_in, map<string, map<string,vector<double> > > &ideal_rna){
    Model_F dbn=read_dbn(filename,1);
    vector<sizet_sizet_vc> gapif_s;
    for(int i=0;i<pdb_in.chains.size();i++){
        vector<int> pairif=pair_info(dbn.chains_F[i]);
        sizet_sizet_vc gapif=gapinfo(fasta_in.chains_F[i],pdb_in.chains[i]);
        gapif_s.push_back(gapif);
        single_pairs(pairif,gapif,fasta_in.chains_F[i],pdb_in.chains[i],ideal_rna);
        vector<int>().swap(pairif);
        sizet_sizet_vc().swap(gapif);
        //cout<<"single_pairs is good!"<<endl;
    }
    return gapif_s;
}

void extrapolate_single_pairs(char* filename, Model_F &fasta_in, ModelUnit &pdb_in){
    Model_F dbn=read_dbn(filename,1);
    for(int i=0;i<pdb_in.chains.size();i++){
        vector<int> pairif=pair_info(dbn.chains_F[i]);
        extrapolate_single_pairs(pairif,fasta_in.chains_F[i],pdb_in.chains[i]);
    }
}


void put_double_pairs(char* filename, Model_F &fasta_in, ModelUnit &pdb_in, vector<sizet_sizet_vc> &gapif_s, vector<pair<double,vector<size_t> > > &bp_vec, map<string, map<string,vector<double> > > &ideal_rna){
    standardize_pdb_order(pdb_in);
    Model_F dbn=read_dbn(filename,1);
    for(int i=0;i<pdb_in.chains.size();i++){
        vector<int> pairif=pair_info(dbn.chains_F[i]);
        //cout<<"double_pairs is bad!"<<endl;
        double_pairs(pairif,gapif_s[i],fasta_in.chains_F[i],pdb_in.chains[i],ideal_rna);
        //cout<<"double_pairs is good!"<<endl;
    }
    standardize_pdb_order(pdb_in);
    MissingRNAatom(pdb_in,ideal_rna,bp_vec,5);
}
void superpose_single_chain(){}

void superpose_chain(){}
#endif