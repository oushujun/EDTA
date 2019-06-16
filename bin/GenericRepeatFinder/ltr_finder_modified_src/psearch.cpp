#include "struct.h"
#include "PBS.h"
#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;

int PBS_minLen=14;
int CHECK_PAIRS=false;

int main(int argc,char* argv[])
{
	if(argc<3)
	{
		cout<<"psearch tRNA.txt seq.txt \n";
		exit(1);
	}
	
	gap_open = 3; 
	gap_ext = 1;
	gap_end = 1;
	score_match = 2;
	score_mismatch = -2;
	set_score_matrix(gap_open, gap_ext, gap_end, score_match, score_mismatch);
	

	
	PBS pbs;
	pbs.LoadSeq(argv[1]);
	ifstream file;
	file.open(argv[2]);
	if(!file.good())
	{
		cerr<<"open "<<argv[2]<<" error\n";
		exit(1);
	}
	string id,str;
	while(file>>id>>str)
	{
		string tmp_name,tmp_pbs;
	    int tmp_begin = 0;
	    int tmp_end = 0;
	    int tmp_match = 0;
	    if(!pbs.Search(str.c_str(), str.length(), tmp_name,
                    tmp_begin, tmp_end, tmp_match, tmp_pbs))
        {
        	cout<<">"<<id<<endl;
        	continue;
        }
        cout<<">"<<id<<" str="<<str<<" name="<<tmp_name<<",begin="<<tmp_begin<<",end="<<tmp_end<<",match_len=";
        cout<<tmp_match<<endl;
        cout<<tmp_pbs<<endl;
    }
		
}
