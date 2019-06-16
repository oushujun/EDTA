/*
* =====================================================================================
*
*       Filename:  PBS.cpp
*
*    Description:  finding PBS
*
*        Version:  1.0
*        Created:  2006年10月22日 23时17分01秒 CST
*       Revision:  none
*       Compiler:  gcc
*
*         Author:   (), 
*        Company:  
*
* =====================================================================================
*/

#include "seq.h"
#include "PBS.h"
#include "struct.h"
#include <stdio.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <map>
#include <string.h>

using namespace std;
void outstr1 (char *str, int pos, int len)
{
    for (int i = pos; i < pos + len && str[i] != '\0'; ++i)
        cout << str[i];
}

char complement_nucleotide(char ch) {
    ch = toupper(ch);

    if ( ch == 'A' )
        return 'T';
    else if ( ch == 'T' )
        return 'A';
    else if ( ch == 'C' )
        return 'G';
    else if ( ch == 'G' )
        return 'C';
    else
        return 'N';
}

PBS::PBS()
{
    init();
}

PBS::PBS(char* filename)
{
    init();
    LoadSeq(filename);
}

void PBS::init()
{
    PBS_len = 18;
    minLen = PBS_minLen;
}

string PBS::LoadSeq(char* filename)
{
    seq_t seq;
    seq.s = NULL;
    seq.m = 0;
    char name[1024];
    char buf[1024];
    int all_count = 0;
    int unic_count = 0;
    string tmp;
    string unic_name;
    map<string, string> unic_seq;
    //cout<<"Load "<<filename<<endl;
    FILE *fp = fopen(filename, "r");

    if (fp == NULL)
    {
        cerr << "open " << filename << " error\n";
        exit (1);
        //return "";
    }
    while ( -1 != read_fasta(fp, &seq, name, 0))
    {
        all_count++;
        tmp = (char*)(seq.s + (seq.l - PBS_len));

        for (int i = 0;i < tmp.length();++i)
            tmp[i] = toupper(tmp[i]);

        //transform(tmp.begin(),tmp.end(),tmp.begin(),toupper);
        //cout<<"len:"<<tmp.length()<<endl;
        //cout<<"seq:"<<name<<"seq len:"<<seq.l<<endl;
        //cout<<"get:"<<tmp<<endl;
        //string tmp1(seq.s);
        //string tmp=tmp1
        int name_len = strlen(name);

        int begin = -1;

        for (int j = 0;j < name_len;++j)
        {
            if (name[j] == '-')
                begin = j;
        }

        unic_name = (char*)(name + (begin + 1));

        if (unic_seq.count(tmp)) //has this seq
        {

            if (unic_seq[tmp].find(unic_name) != string::npos)
                continue;
            else
            {
                unic_seq[tmp] += "|";
                unic_seq[tmp] += unic_name;
            }

        }
        else
            unic_seq[tmp] = unic_name;
    }

    free(seq.s);

    for (map<string, string>::iterator it = unic_seq.begin(); it != unic_seq.end(); ++it)
    {
        unic_count++;
        tmp = it->first;
        _tRNA.push_back(tmp);
        //cout<<">"<<it->second<<endl;
        //cout<<tmp<<endl;
        //reverse(_tRNA.rbegin()->begin(), _tRNA.rbegin()->end());
        reverse(tmp.begin(), tmp.end());
        transform(tmp.begin(), tmp.end(), tmp.begin(),
                  complement_nucleotide); //get complement seq
        tRNA.push_back( tmp );
        //cout<<tmp<<endl;
        unic_name = it->second;
        tRNA_name.push_back(unic_name);
        _tRNA_name.push_back((string)"-" + unic_name);
    }

    fclose(fp);
    sprintf(buf, "(unic:%d/all:%d) ", unic_count, all_count);
    return string(buf);
}

int PBS::Search(const char* str, int len, string& name,
                int& begin, int& end, int& match_len, string& outstr)
{
    if (len < PBS_len)
        return 0;

    int path_len = 0;

    int max = -100000;

    int mLen = 0;

    int index = -1;

    string n;

    int b = 0, e = 0, m = 0;

    for (int i = 0;i < tRNA.size();++i)
    {
        //cout<<i<<" begin match:"<<tRNA[i]<<endl;

        AlnAln *aln;
        aln = aln_stdaln_aux(str, tRNA[i].c_str(), &aln_param_nt2nt, 0,
                             len, PBS_len);

        //cout<<"end  match"<<endl;
        //    cout<<"SEQ:";outstr1(str,0,len);cout<<endl;
        //    cout<<"RNA:"<<tRNA[i]<<endl;

        if (aln->path_len > 0 && aln->score > max)
        {
            max = aln->score;
            index = i;
            n = tRNA_name[i];
            b = aln->start1 - 1;
            e = aln->path_len;
            m = 0;

            for (int j = 0;j < aln->path_len;++j)
                if (aln->outm[j] == '|')
                    ++m;

            outstr = "";

            outstr += aln->out2;

            outstr += "\n";

            outstr += aln->outm;

            outstr += "\n";

            outstr += aln->out1;

            //cout<<i<<" PBS: "<<name<<" match_len: "<<m<<" score:"<<max<<endl;
            //cout<<"     "<<aln->out1<<endl;
            //cout<<"     "<<aln->outm<<endl;
            //cout<<"     "<<aln->out2<<endl;

        }

        aln_free_AlnAln(aln);
    }
    if (m >= minLen)
    {
        name = n;
        begin = b;
        end = b + e - 1;
        match_len = m;
        return max;
    }
    else
        return 0;
}

int PBS::rSearch(const char* str, int len, string& name,
                 int& begin, int& end, int& match_len, string& outstr)
{
    if (len < PBS_len)
        return 0;

    int path_len = 0;

    int max = -100000;

    int mLen = 0;

    int index = -1;

    string n;

    int b = 0, e = 0, m = 0;


    for (int i = 0;i < _tRNA.size();++i)
    {

        AlnAln *aln;
        aln = aln_stdaln_aux(str, _tRNA[i].c_str(), &aln_param_nt2nt, 0,
                             len, PBS_len);

        if (aln->path_len > 0 && aln->score > max)
        {
            max = aln->score;
            index = i;
            n = _tRNA_name[i];
            b = aln->start1 - 1;
            e = aln->path_len;
            m = 0;

            for (int j = 0;j < aln->path_len;++j)
                if (aln->outm[j] == '|')
                    ++m;

            outstr = "";

            outstr += aln->out2;

            outstr += "\n";

            outstr += aln->outm;

            outstr += "\n";

            outstr += aln->out1;

            //cout<<"PBS: "<<name<<" match_len: "<<match_len<<endl;
            //cout<<"     "<<out1<<endl;
            //cout<<"     "<<out2<<endl;
        }

        aln_free_AlnAln(aln);
    }
    if (m >= minLen)
    {
        name = n;
        begin = b;
        end = b + e - 1;
        match_len = m;
        return max;
    }
    else
        return 0;
}


//////////////////////////////////////
PPT::PPT()
{
    init();
}
PPT::PPT(int w, int d)
{
    window_size = w; //reset window size
    delta = d;

}
void PPT::init()
{
    window_size = 17;
    delta = 3;
}


int PPT::Search(char *string, int len, int& begin, int& match_count)
{
    match_count = 0;

    if (len < window_size) //region too short
        return 0;

    int b = 0;

    int target_count = 0;

    //first window  // need to rewrite, go backword
    for (int i = 0;i < window_size;++i)
    {
        if (string[i] == 'A' || string[i] == 'G')
            ++target_count;
    }

    int max = target_count;

    for (int i = 1;i <= len - window_size;++i)
    {
        if (string[i - 1] == 'A' ||
                string[i - 1] == 'G')
            --target_count;

        if (string[i + window_size - 1] == 'A' ||
                string[i + window_size - 1] == 'G')
            ++target_count;

        if (max <= target_count) //right max
        {
            max = target_count;
            b = i;
        }

    }

    if (window_size - max <= delta )
    {
        begin = b;
        match_count = max;
        return 1;
    }

    return 0;
}

int PPT::rSearch(char *string, int len, int& begin, int& match_count)
{
    match_count = 0;

    if (len < window_size) //region too short
        return 0;

    int b = 0;

    int target_count = 0;

    //first window
    for (int i = 0;i < window_size;++i)
    {
        if (string[i] == 'T' || string[i] == 'C')
            ++target_count;
    }

    int max = target_count;

    for (int i = 1;i <= len - window_size;++i)
    {
        if (string[i - 1] == 'T' ||
                string[i - 1] == 'C')
            --target_count;

        if (string[i + window_size - 1] == 'T' ||
                string[i + window_size - 1] == 'C')
            ++target_count;

        if (max < target_count) //left max
        {
            max = target_count;
            b = i;
        }

    }

    if (window_size - max <= delta)
    {
        begin = b;
        match_count = max;
        return 1;
    }

    return 0;
}


