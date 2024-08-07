#include "PDBParser.h"
#include "cssr_struct.h"

const char * docstring=""
"CSSR input.fasta input.pdb output.dbn\n"
"\n"
"Input\n"
"    input.fasta - full length sequence\n"
"    input.pdb   - input structure. residue index must be consistent with input.fasta\n"
"\n"
"Output\n"
"    output.dbn  - dot bracket format secondary structure\n"
;

int main( int argc, char* argv[] ) 
{
    if (argc!=4)
    {
        cerr<<docstring<<endl;
        return 0;
    }
    string fastafile = argv[1];
    string input     = argv[2];
    string ctFile    = argv[3];


    // parse 3D structure
    int atomic_detail=2;
    int allowX       =1; // only allow ATOM and MSE
    ios_base::sync_with_stdio(false);
    ModelUnit pdb_entry=read_pdb_structure(input.c_str(),atomic_detail,allowX);
    
    // parse sequence into upper case
    string sequence="";
    ifstream fin;
    string line;
    if (fastafile!="-") fin.open(fastafile.c_str());
    while ((fastafile=="-")?cin.good():fin.good())
    {
        if (fastafile=="-") getline(cin,line);
        else getline(fin, line);
        if (line.size() && line[0]!='>') sequence+=line;
    }
    fin.close();
    line.clear();

    // assign SS solely by structure
    vector<string> res_str_vec;
    vector<pair<float,vector<string> > > bp_vec;
    bool interchain=0;
    cssr(pdb_entry, res_str_vec, bp_vec, interchain);
    size_t bp;
    vector<size_t> filtered_bp_vec;
    filter_bp(bp_vec, filtered_bp_vec);

    ofstream fp;
    bool usestdout=(ctFile=="-");
    if (!usestdout) fp.open(ctFile.c_str());
    vector<char> dot_bracket;
    cssr_bp2dot(res_str_vec, filtered_bp_vec, bp_vec, dot_bracket);
    
    // map SS to sequence
    vector<char> dot_bracket_fill(sequence.size(), '.');
    size_t c,r;
    int resi;
    size_t L=0;
    for (c=0; c<pdb_entry.chains.size(); c++)
    {
        for (r=0; r<pdb_entry.chains[c].residues.size(); r++)
        {
            resi=pdb_entry.chains[c].residues[r].resi;
            if (resi<=0 || resi>sequence.size()) continue;
            dot_bracket_fill[resi-1]=dot_bracket[r+L];
        }
        L+=pdb_entry.chains[c].residues.size();
    }


    if (usestdout)
    {
        cout<<sequence<<endl;
        for (size_t r=0;r<dot_bracket_fill.size();r++) cout<<dot_bracket_fill[r];
        cout<<endl;
    }
    else
    {
        fp<<sequence<<endl;
        for (size_t r=0;r<dot_bracket_fill.size();r++) fp<<dot_bracket_fill[r];
        fp<<endl;
    }
    dot_bracket.clear();
    if (!usestdout) fp.close();

	/* clean up */
    vector<size_t>().swap(filtered_bp_vec);
    vector<string>().swap(res_str_vec);
    vector<pair<float,vector<string> > >().swap(bp_vec);
    vector<ChainUnit>().swap(pdb_entry.chains);
    return 1;
}

