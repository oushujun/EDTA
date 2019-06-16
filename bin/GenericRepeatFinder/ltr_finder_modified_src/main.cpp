#include "struct.h"
#include "PairsFilter.h"
#include "seq.h"
#include "LinearSuffixSort.h"
#include "PBS.h"
#include "CallPsScan.h"
#include "Timing.h"
#include "regex.h"

#include <stdio.h>
#include <cstring>
#include <vector>
#include <algorithm>
#include <fstream>
#include <unistd.h>
#include <set>
#include <sstream> 
#include <map>

using namespace std;

//global variable define
int Dmin = 1000;
int Dmax = 20000;
int Lmin = 100;
int Lmax = 3500;
int Lex = 20; //exact match pair len
int MaxGap = 50; //maxgap between pairs
float JoinThreshold = 0.7;
float SplitThreshold = 0.9;
int TSRmin = 4;
int TSRmax = 6;
int PBS_region = 100;
int PPT_region = 100;
int PPT_window = 15;
int PBS_minLen = 14;
int max_sub_rt_gap = 2;
int min_sub_rt_count = 4;
int CHECK_PAIRS = 0;
float minOutScore = 6.0;
int wrought = 0;
float LowerSharpness = 0.40;
float HigherSharpness = 0.40;
float minMatchSim = 0;
int outAlignLen = 40;
char *fig_file = NULL;
bool edge_signal = false;
bool showPairNum = false;
bool checkCentriole = false;
char transDNA[128];
char transAA[128];
char* Filter = NULL;
//5'TG 5'CA 3'TG 3'CA TSR PBS PPT RTdomain Core C-Term
//  1   1     1   1    1   1   1     1      1     1
//  1111111111B
int SignalNum = 10;
int FModule = (1 << (SignalNum + 1)) - 1; //module
int M5TG = 1 << 10; //10000000000b;
int M5CA = 1 << 9; //01000000000b;
int M3TG = 1 << 8; //00100000000b;
int M3CA = 1 << 7; //00010000000b;
int MTSR = 1 << 6; //00001000000b;
int MPBS = 1 << 5; //00000100000b;
int MPPT = 1 << 4; //00000010000b;
int MRT = 1 << 3; //00000001000b;
int MCORE = 1 << 2;
int MCT = 1 << 1;
int MRH = 1;

string matrix2domain(const string& m) {
    //    if(m.find("PS50878") != string::npos)
    //    {
    //        return "RT Domain";
    //    }
    //    else

    if (m.find("PS50994") != string::npos) {
        return "IN (core)";
    } else if (m.find("PS51027") != string::npos) {
        return "IN (c-term)";
    } else if (m.find("PS50879") != string::npos) {
        return "RNase H";
    } else {
        return m;
    }
}

int matrix2modul(string& m) {
    if (m.find("RT") != string::npos) {
        return MRT;
    } else if (m.find("PS50994") != string::npos) {
        return MCORE;
    } else if (m.find("PS51027") != string::npos) {
        return MCT;
    } else if (m.find("PS50879") != string::npos) {
        return MRH;
    } else {
        return 0;
    }

}

void CountScore(stick& st) {


    string PBS_name, PBS_str, f_out, r_out;
    int PBS_begin, PBS_end, PPT_begin, PPT_count;
    char buf[10240]; //2014-07-30, increase buf size to prvent segmentation fault, by Zhao

    //try to output foward strand
    int f_score = 0;
    int f_dscore = 0;
    int f_status = 0;
    int f_pp_match = 0; //PBS and PPT match base count
    int a_begin = st.end1 + 1;
    int a_end = st.pos2 - 1;

    if (st.pp.PBS) {
        f_score++;
        f_status = f_status | MPBS;
        //            sprintf(buf,"PBS   : [%s] %d - %d Len: %d\n",st.pp.PBS_name.c_str(),
        //                    st.pp.PBS_begin+1,st.pp.PBS_end+1,
        //                    st.pp.PBS_match);
        sprintf(buf, "PBS   : [%d/%d] %d - %d (%s)\n", st.pp.PBS_match,
                st.pp.PBS_end - st.pp.PBS_begin + 1,
                st.pp.PBS_begin + 1, st.pp.PBS_end + 1, st.pp.PBS_name.c_str());
        f_out += buf;
        a_begin = st.pp.PBS_end + 1;

        for (int i = 0; i < st.pp.pbs_str.length(); ++i) {
            if (st.pp.pbs_str[i] == '|')
                ++f_pp_match;
        }

    }

    if (st.pp.PPT) {
        f_score++;
        f_status = f_status | MPPT;
        sprintf(buf, "PPT   : [%d/%d] %d - %d\n", st.pp.PPT_count, PPT_window,
                st.pp.PPT_begin + 1, st.pp.PPT_begin + PPT_window);
        f_out += buf;
        a_end = st.pp.PPT_begin - 1;
        f_pp_match += st.pp.PPT_count;
    }
    for (int j = 0; j < st.motif.size(); ++j) {
        if (st.motif[j].forward) {
            int p_begin = st.motif[j].p_begin[0]; //ORF: from * to *
            int p_end = st.motif[j].p_end[st.motif[j].p_end.size() - 1];
            //                int p_begin=st.motif[j].begin;
            //                int p_end=st.motif[j].end;
            //                for(int k=0;k<st.motif[j].p_begin.size();++k)
            //                {
            //                    if(st.motif[j].p_begin[k]>=a_begin &&
            //                            st.motif[j].p_begin[k]<=a_end)
            //                    {
            //                        p_begin=st.motif[j].p_begin[k];
            //                        break;
            //                    }
            //                }
            //                for(int k=(int)st.motif[j].p_end.size()-1;k>=0;--k)
            //                {
            //                    if(st.motif[j].p_end[k]>=a_begin &&
            //                            st.motif[j].p_end[k]<=a_end)
            //                    {
            //                        p_end=st.motif[j].p_end[k];
            //                        break;
            //                    }
            //                }
            sprintf(buf, "Domain: %d - %d [possible ORF:%d-%d, (%s)]\n",
                    st.motif[j].begin + 1, st.motif[j].end + 1, p_begin + 1, p_end + 1,
                    matrix2domain(st.motif[j].name).c_str());
            f_out += buf;
            f_score++;
            f_dscore++;
            f_status = f_status | matrix2modul(st.motif[j].name);
        }

    }

    //try to output reverse strand
    int r_score = 0;
    int r_dscore = 0;
    int r_status = 0;
    int r_pp_match = 0;
    a_begin = st.end1 + 1;
    a_end = st.pos2 - 1;

    if (st.pp._PBS) {
        r_score++;
        r_status = r_status | MPBS;
        //            sprintf(buf,"PBS   : [%s] %d - %d Len: %d\n",st.pp._PBS_name.c_str(),
        //                    st.pp._PBS_begin+1,st.pp._PBS_end+1,
        //                    st.pp._PBS_match);
        sprintf(buf, "PBS   : [%d/%d] %d - %d (%s)\n", st.pp._PBS_match,
                st.pp._PBS_end - st.pp._PBS_begin + 1,
                st.pp._PBS_begin + 1, st.pp._PBS_end + 1, st.pp._PBS_name.c_str());
        r_out += buf;
        a_end = st.pp._PBS_begin - 1;

        for (int i = 0; i < st.pp._pbs_str.length(); ++i) {
            if (st.pp._pbs_str[i] == '|')
                ++r_pp_match;
        }

    }

    if (st.pp._PPT) {
        r_score++;
        r_status = r_status | MPPT;
        sprintf(buf, "PPT   : [%d/%d] %d - %d\n", st.pp._PPT_count, PPT_window,
                st.pp._PPT_begin + 1, st.pp._PPT_begin + PPT_window);
        r_out += buf;
        a_begin = st.pp._PPT_begin + st.pp._PPT_count;
        r_pp_match += st.pp._PPT_count;
    }
    for (int j = 0; j < st.motif.size(); ++j) {
        if (!st.motif[j].forward) {
            int p_begin = st.motif[j].p_begin[0];
            int p_end = st.motif[j].p_end[st.motif[j].p_end.size() - 1];
            //                int p_begin=st.motif[j].begin;
            //                int p_end=st.motif[j].end;
            //                for(int k=0;k<st.motif[j].p_begin.size();++k)
            //                {
            //                    if(st.motif[j].p_begin[k]>=a_begin &&
            //                            st.motif[j].p_begin[k]<=a_end)
            //                    {
            //                        p_begin=st.motif[j].p_begin[k];
            //                        break;
            //                    }
            //                }
            //                for(int k=(int)st.motif[j].p_end.size()-1;k>=0;--k)
            //                {
            //                    if(st.motif[j].p_end[k]>=a_begin &&
            //                            st.motif[j].p_end[k]<=a_end)
            //                    {
            //                        p_end=st.motif[j].p_end[k];
            //                        break;
            //                    }
            //                }
            sprintf(buf, "Domain: %d - %d [possible ORF:%d-%d, (%s)]\n",
                    st.motif[j].begin + 1, st.motif[j].end + 1, p_begin + 1, p_end + 1,
                    matrix2domain(st.motif[j].name).c_str());
            r_out += buf;
            r_score++;
            r_dscore++;
            r_status = r_status | matrix2modul(st.motif[j].name);
        }

    }
    /* //I think this program should not output TWO STRAND LTR
    if( (f_score>=3 && r_score>=3) ||
            (f_score==r_score && f_score>0) )
    {
        st.score+=f_score;
        st.status = st.status | f_status;
        st.strand_str+="BOTH STRANDS ARE POSSIBLE:\n";
        st.strand_str+="Strand + :\n";
        st.strand_str+=f_out;
        st.strand_str+="Strand - :\n";
        st.strand_str+=r_out;
        st.strand=3;//two
    }
    else */

    if (f_score > r_score ||
            ((f_score == r_score) && (f_pp_match >= r_pp_match) && (f_score > 0))) {
        st.score += f_score;
        st.domain_score = f_dscore;
        st.status = st.status | f_status;
        st.strand_str += "Strand + :\n";
        st.strand_str += f_out;
        st.strand = 1; //+
    } else if (f_score < r_score ||
            ((f_score == r_score) && (f_pp_match < r_pp_match) && (r_score > 0))) {
        st.score += r_score;
        st.domain_score = r_dscore;
        st.status = st.status | r_status;
        st.strand_str += "Strand - :\n";
        st.strand_str += r_out;
        st.strand = 2; //-
    } else {
        st.strand_str += "NO EVIDENCE FOR STRAND\n";
        st.strand = 0; //no
    }

}

string dec2bin(int dec) {
    string ret;

    for (int i = SignalNum; i >= 0; --i) {
        if (dec & (1 << i))
            ret += "1";
        else
            ret += "0";
    }

    return ret;
}

int InitFilter(char* f) {
    if (f == 0)
        return FModule;

    int ret = f[0] - '0';

    for (int i = 1; i < strlen(f); ++i) {
        ret = ret << 1;
        ret += f[i] - '0';
    }

    return ~ret & FModule;
}

void TableOutputMotif(vector<MOTIF>& motif, string name, bool forward) {
    bool has = false;

    for (int j = 0; j < motif.size(); ++j) {
        if (motif[j].forward == forward) {
            if (motif[j].name == name) {
                cout << motif[j].begin + 1 << "-" << motif[j].end + 1 << "\t";
                has = true;
                break;
            }

        }

    }

    if (!has)
        cout << "N-N\t";
}

bool FigOutputMotif(vector<MOTIF>& motif, string name, bool forward, ofstream& fig, string id) {
    bool has = false;

    for (int j = 0; j < motif.size(); ++j) {
        if (motif[j].forward == forward) {
            if (motif[j].name == name) {
                fig << motif[j].begin + 1 << "," << motif[j].end + 1 << "," << id << ";";
                has = true;
                //break;
            }

        }

    }
    return has;
}

int OutPutResult(char* name, char* myString, int stringLength, vector < stick >& sticks, ofstream& fig) {
    int result_count = 0;
    FourPos fp;

    set
            < FourPos > used_ltr;

    int filter = InitFilter(Filter);

    if (fig.good()) {
        fig << ">" << name << "\n";
    }

    cout << ">Sequence: " << name << " Len:" << stringLength << endl;

    if (wrought == 2)
        cout << "index\tSeqID\tLocation\tLTR len\tInserted element len\tTSR\tPBS\tPPT\tRT\tIN (core)\tIN (c-term)\tRH\tStrand\tScore\tSharpness\tSimilarity\n";

    // I think we should not del overlap by hand, let judge it later
    // And, if we del these kind of overlap, we should check each two of them,
    // not only the neighbored.
    //    //del overlap
    //    for(int i=0;i<(int)sticks.size()-1 && Sensitive<0.001;++i)
    //    {
    //        cerr<<"check :"<<sticks[i].pos1+1<<"-"<<sticks[i].end2+1<<"  "
    //            <<sticks[i+1].pos1+1<<"-"<<sticks[i+1].end2+1<<endl;
    //        if(sticks[i].end2<sticks[i+1].pos1)
    //            continue;//no overlap
    //        int l_len=0;
    //        int r_len=0;
    //        int m_len=0;
    //        if(sticks[i].end1>sticks[i+1].pos1 &&
    //                sticks[i].pos1<sticks[i+1].end1 )
    //        {
    //            int max = sticks[i].end1>sticks[i+1].end1?sticks[i].end1:sticks[i+1].end1;
    //            int min = sticks[i].pos1<sticks[i+1].pos1?sticks[i].pos1:sticks[i+1].pos1;
    //            l_len=max-min;
    //            int min_ltr = sticks[i].len1<sticks[i+1].len1?sticks[i].len1:sticks[i+1].len1;
    //            if( (float)l_len/min_ltr < 0.9)
    //                l_len=0;//this should not regard as overlap
    //        }
    //        if(sticks[i].end2>sticks[i+1].pos2 &&
    //                sticks[i].pos2<sticks[i+1].end2)
    //        {
    //            int max = sticks[i].end2>sticks[i+1].end2?sticks[i].end2:sticks[i+1].end2;
    //            int min = sticks[i].pos2<sticks[i+1].pos2?sticks[i].pos2:sticks[i+1].pos2;
    //            r_len=max-min;
    //            int min_ltr = sticks[i].len2<sticks[i+1].len2?sticks[i].len2:sticks[i+1].len2;
    //            if( (float)r_len/min_ltr < 0.9)
    //                r_len=0;//this should not regard as overlap
    //        }
    //        if(sticks[i].end2>sticks[i+1].pos1 &&
    //                sticks[i].pos2<sticks[i+1].end1)
    //        {
    //            int max = sticks[i].end2>sticks[i+1].end1?sticks[i].end2:sticks[i+1].end1;
    //            int min = sticks[i].pos2<sticks[i+1].pos1?sticks[i].pos2:sticks[i+1].pos1;
    //            m_len=max-min;
    //            int min_ltr = sticks[i].len2<sticks[i+1].len1?sticks[i].len2:sticks[i+1].len1;
    //            if( (float)m_len/min_ltr < 0.9)
    //                m_len=0;//this should not regard as overlap
    //        }
    //        if(l_len>0 && r_len>0)//both LTR overlap
    //        {
    //            if(sticks[i].score < sticks[i+1].score)
    //                sticks[i].score = 0;//del the lower score, if they are same, no one deleted
    //            else if(sticks[i].score > sticks[i+1].score)
    //                sticks[i+1].score = 0;
    //        }
    //        else if(l_len>0 || r_len>0)//only one side same
    //        {
    //            cerr<<" one side equal, domain:"<<sticks[i].motif.size()<<" "<<sticks[i+1].motif.size()<<endl;
    //            //should del too more aa one, or no aa one
    //            if(sticks[i].motif.size()>0 && sticks[i].motif == sticks[i+1].motif)//domain same
    //            {
    //                cerr<<"motif same\n";
    //                if(sticks[i].score < sticks[i+1].score)
    //                    sticks[i].score = 0;//del the lower score, if they are same, no one deleted
    //                else if(sticks[i].score > sticks[i+1].score)
    //                    sticks[i+1].score = 0;
    //            }
    //
    //        }
    //        else if(m_len>0)//first 3'LTR overlap with second 5'LTR
    //        {
    //            if(sticks[i].domain_score>0 && sticks[i+1].domain_score<=0)//one has no domain
    //                sticks[i+1].score=0;//del
    //            else if(sticks[i].domain_score<=0 && sticks[i+1].domain_score>0)//one has no domain
    //                sticks[i].score=0;//del
    //        }
    //    }

    for (int i = 0; i < sticks.size(); ++i) {
        stick st = sticks[i];

        if ((!CHECK_PAIRS && st.score < minOutScore) || (CHECK_PAIRS && st.score < 0))
            continue; //

        fp.a = st.pos1;

        fp.b = st.end1;

        fp.c = st.pos2;

        fp.d = st.end2;

        if (used_ltr.count(fp))
            continue;
        else
            used_ltr.insert(fp);

        //recheck condicitons
        if (st.len1 < Lmin || st.len1 > Lmax ||
                st.len2 < Lmin || st.len2 > Lmax ||
                (st.pos2 - st.end1) < Dmin || (st.pos2 - st.end1) > Dmax) {
            //cout<<"filter, this one:"<<st.pos1<<"-"<<st.end2<<endl;
            continue;
        }

        //check sharpness
        if (st.sharpness5 < LowerSharpness || st.sharpness3 < LowerSharpness)
            continue;

        if (st.sharpness5 < HigherSharpness && st.sharpness3 < HigherSharpness)
            continue;

        //if(st.match_score<0.83)
        //    continue;
        //check Filter
        if (FModule != (st.status | filter)) {
            //cout<<st.status<<" filter:"<<filter<<endl;
            //cout<<"ans: "<<(st.status|filter)<<" M:"<<FModule<<endl;
            continue;
        }

        //        if (st.match_score < minMatchSim)
        //            continue;

        if (edge_signal) //check edge signal
        {

            if (st.strand == 1) {
                if (st.pp.PBS + st.pp.PPT + (st.tsr_len > 0) < 2)
                    continue;
            } else if (st.strand == 2) {
                if (st.pp._PBS + st.pp._PPT + (st.tsr_len > 0) < 2)
                    continue;
            } else
                continue; //no strand? not output
        }

        result_count++;

        if (fig.good()) {
            fig << st.pos1 + 1 << "," << st.end1 + 1 << ";" << st.pos2 + 1 << "," << st.end2 + 1 << "\t";
            fig << setprecision(3);
            fig << "[" << result_count << "]LTR(" << " +-?"[st.strand] << ") score:" << st.score << "(" << st.match_score << ")\t";
            FigOutputMotif(st.motif, "RT", st.strand == 1, fig, "RT"); //RT
            FigOutputMotif(st.motif, "PS50994", st.strand == 1, fig, "IN(core)"); //IN Core
            FigOutputMotif(st.motif, "PS51027", st.strand == 1, fig, "IN(c-term)"); //IN C term
            FigOutputMotif(st.motif, "PS50879", st.strand == 1, fig, "RH"); //RNaseH

            if (st.strand == 1 && st.pp.PBS)
                fig << st.pp.PBS_begin + 1 << "," << st.pp.PBS_end + 1 << ",PBS;";
            else if (st.strand == 2 && st.pp._PBS)
                fig << st.pp._PBS_begin + 1 << "," << st.pp._PBS_end + 1 << ",PBS;";

            if (st.strand == 1 && st.pp.PPT)
                fig << st.pp.PPT_begin + 1 << "," << st.pp.PPT_begin + PPT_window << ",PPT;";
            else if (st.strand == 2 && st.pp._PPT)
                fig << st.pp._PPT_begin + 1 << "," << st.pp._PPT_begin + PPT_window << ",PPT;";

            /* output flank bases
            int begin=st.pos1-5;
            int end=st.pos1-1;
            if(begin<0) begin=0;
            if(end>strlen(myString)) end=strlen(myString)-1;
            fig<<begin+1<<","<<end+1<<",";
            for(int k=begin;k<=end;++k)
                fig<<myString[k];
            fig<<"|"<<myString[st.pos1]<<myString[st.pos1+1];
            fig<<";";
            begin=st.end2+1;
            end=st.end2+5;
            if(begin<0) begin=0;
            if(end>strlen(myString)) end=strlen(myString)-1;
            fig<<begin+1<<","<<end+1<<",";
            fig<<myString[st.end2-1]<<myString[st.end2]<<"|";
            for(int k=begin;k<=end;++k)
                fig<<myString[k];
            fig<<";";
             */
            //output TSR
            if (st.tsr_len) {
                fig << st.pos1 - st.tsr_len + 1 << "," << st.pos1 - 1 + 1 << ",TSR;";
                fig << st.end2 + 1 + 1 << "," << st.end2 + st.tsr_len + 1 << ",TSR;";
            }

            fig << endl;
        }
        if (wrought == 2) //table output
        {
            string TSR = "N";

            if (st.tsr_len) {
                TSR = "";

                for (int i = st.tsr_pos1; i < st.tsr_pos1 + st.tsr_len; ++i)
                    TSR += myString[i];
            }

            char buf[1024];
            //index,name,LTR region,LTR len,TG,CA,TSR,Score
            sprintf(buf, "[%2d]\t%s\t%d-%d\t%d,%d\t%d\t%s\t",
                    result_count, name, st.pos1 + 1, st.end2 + 1, st.len1, st.len2, st.end2 - st.pos1 + 1,
                    TSR.c_str());
            cout << buf;

            if (st.strand == 1) //forward strand
            {

                if (st.pp.PBS)
                    cout << st.pp.PBS_begin + 1 << "-" << st.pp.PBS_end + 1 << "\t";
                else
                    cout << "N-N\t";

                if (st.pp.PPT)
                    cout << st.pp.PPT_begin + 1 << "-" << st.pp.PPT_begin + PPT_window << "\t";
                else
                    cout << "N-N\t";
            } else if (st.strand == 2) {
                if (st.pp._PBS)
                    cout << st.pp._PBS_begin + 1 << "-" << st.pp._PBS_end + 1 << "\t";
                else
                    cout << "N-N\t";

                if (st.pp._PPT)
                    cout << st.pp._PPT_begin + 1 << "-" << st.pp._PPT_begin + PPT_window << "\t";
                else
                    cout << "N-N\t";
            } else {
                cout << "N-N\tN-N\t";
            }

            TableOutputMotif(st.motif, "RT", st.strand == 1); //RT
            TableOutputMotif(st.motif, "PS50994", st.strand == 1); //IN Core
            TableOutputMotif(st.motif, "PS51027", st.strand == 1); //IN C term
            TableOutputMotif(st.motif, "PS50879", st.strand == 1); //RNaseH

            if (st.strand == 1)
                cout << "+\t";
            else if (st.strand == 2)
                cout << "-\t";
            else
                cout << "?\t";

            cout << setprecision(3);

            cout << st.score << "\t" << st.sharpness5 << "," << st.sharpness3;

            cout << setprecision(3);

            cout << "\t" << st.match_score;

            cout << endl;

            continue;
        }

        cout << "[" << result_count << "] " << name << " Len:" << stringLength << endl;
        cout << "Location : " << st.pos1 + 1 << " - " << st.end2 + 1 << " Len: " << st.end2 - st.pos1 + 1 << " Strand:";

        if (st.strand == 1)
            cout << "+" << endl;
        else if (st.strand == 2)
            cout << "-" << endl;
        else if (st.strand == 3)
            cout << "+/-" << endl;
        else
            cout << "?" << endl;

        cout << setprecision(3);

        cout << "Score    : " << st.score << " [LTR region similarity:" << st.match_score << "]" << endl;

        cout << "Status   : " << dec2bin(st.status) << endl;

        cout << "5'-LTR   : " << st.pos1 + 1 << " - " << st.end1 + 1 << " Len: " << st.len1 << endl;

        cout << "3'-LTR   : " << st.pos2 + 1 << " - " << st.end2 + 1 << " Len: " << st.len2 << endl;

        if (st.tg_pos1 >= 0) {
            cout << "5'-TG    : " << myString[st.tg_pos1] << myString[st.tg_pos1 + 1] << " , "
                    << myString[st.tg_pos2] << myString[st.tg_pos2 + 1] << endl;
        } else
            cout << "5'-TG    : NOT FOUND\n";

        if (st.ca_pos2 >= 0) {
            cout << "3'-CA    : " << myString[st.ca_pos1] << myString[st.ca_pos1 + 1] << " , "
                    << myString[st.ca_pos2] << myString[st.ca_pos2 + 1] << endl;
        } else
            cout << "3'-CA    : NOT FOUND\n";

        if (st.tsr_len) {
            cout << "TSR      : " << st.tsr_pos1 + 1 << " - "
                    << st.tsr_pos1 + st.tsr_len << " , "
                    << st.tsr_pos2 + 1 << " - "
                    << st.tsr_pos2 + st.tsr_len //<<" Len:"<<st.tsr_len<<endl;
                    << " [";
            outstr(myString, st.tsr_pos1, st.tsr_len);
            cout << "]" << endl;
        } else {
            cout << "TSR      : NOT FOUND" << endl;
        }

        cout << setprecision(3);
        cout << "Sharpness: " << st.sharpness5 << "," << st.sharpness3 << endl;
        cout << st.strand_str << endl;

        if (wrought)
            continue;

        //output pairs
        cout << "Details of exact match pairs:\n";

        for (int j = 0; j < st.candi.size(); ++j) {
            if (j > 0)
                cout << " (" << st.candi[j].pos2 - (st.candi[j - 1].pos2 + st.candi[j - 1].len) << ") ";

            cout << st.candi[j].pos2 + 1 << "-" << st.candi[j].pos2 + st.candi[j].len << "[" << st.candi[j].len << "]";
        }

        cout << endl;

        for (int j = 0; j < st.candi.size(); ++j) {
            if (j > 0)
                cout << " (" << st.candi[j].pos1 - (st.candi[j - 1].pos1 + st.candi[j - 1].len) << ") ";

            cout << st.candi[j].pos1 + 1 << "-" << st.candi[j].pos1 + st.candi[j].len << "[" << st.candi[j].len << "]";
        }

        cout << endl << endl;

        //output detail
        cout << "Details of the LTR alignment(5'-end):" << endl;
        cout << st.LTR5;
        cout << "Details of the LTR alignment(3'-end):" << endl;
        cout << st.LTR3;

        if (st.strand & 1) //1 or 3
        {

            if (st.pp.PBS) {
                cout << "Details of the PBS alignment(+):" << endl;
                cout << "tRNA type: " << st.pp.PBS_name << endl;
                cout << st.pp.pbs_str << endl;
                cout << "|" << st.pp.PBS_begin + 1 << endl << endl;
            }
            if (st.pp.PPT) {
                cout << "Details of PPT(+):\n";
                outstr(myString, st.pp.PPT_begin, PPT_window);
                cout << endl;
                cout << "|" << st.pp.PPT_begin + 1 << endl << endl;
            }

        }

        if (st.strand & 2) //1 or 3
        {

            if (st.pp._PBS) {
                cout << "Details of the PBS alignment(-):" << endl;
                cout << "tRNA type: " << st.pp._PBS_name << endl;
                cout << st.pp._pbs_str << endl;
                cout << "|" << st.pp._PBS_begin + 1 << endl << endl;
            }
            if (st.pp._PPT) {
                cout << "Details of PPT(-):\n";
                outstr(myString, st.pp._PPT_begin, PPT_window);
                cout << endl;
                cout << "|" << st.pp._PPT_begin + 1 << endl << endl;
            }

        }
        //cout<<endl;
    }

    if (result_count == 0)
        cout << "No LTR Retrotransposons Found\n";
    
    return result_count;
}

int main(int argc, char *argv[]) {
    char *myString;
    int stringLength;
    int i;
    string namePattern = "";
    char *version = "1.06";

    bool display_usage = false;
    char c;
    CPSSCAN ps_scan;
    PBS pbs;
    PPT ppt(PPT_window, 4);
    char *tRNA_file = NULL;
    char *ps_dir = NULL;
    bool ishtml = false;

    struct re_pattern_buffer id_filter;
    id_filter.allocated = 0;
    id_filter.buffer = 0;
    id_filter.fastmap = 0;
    id_filter.translate = 0;

    Timing timehere;
    Timing timeall;
    timeall.markbeg();

    while ((c = getopt(argc, argv, "ho:t:e:m:u:D:d:L:l:p:s:cw:S:a:P:g:F:B:b:J:j:O:r:M:f:xEiCG:T:")) > 0) {
        switch (c) {

            case 'h':
                display_usage = true;
                break;

            case 'o':
                gap_open = atoi(optarg);
                break;

            case 't':
                gap_ext = atoi(optarg);
                break;

            case 'e':
                gap_end = atoi(optarg);
                break;

            case 'm':
                score_match = atoi(optarg);
                break;

            case 'u':
                score_mismatch = atoi(optarg);
                break;

            case 'D':
                Dmax = atoi(optarg);
                break;

            case 'd':
                Dmin = atoi(optarg);
                break;

            case 'L':
                Lmax = atoi(optarg);
                break;

            case 'l':
                Lmin = atoi(optarg);
                break;

            case 'p':
                Lex = atoi(optarg);
                break;

            case 'c':
                CHECK_PAIRS = 1;
                break;

            case 'x':
                ishtml = true;
                break;

            case 'E':
                edge_signal = true;
                break;

            case 'i':
                showPairNum = true;
                break;

            case 'w':
                wrought = atoi(optarg);
                break;

            case 'P':
                namePattern = optarg;
                break;

            case 'F':
                Filter = optarg;
                break;

            case 'g':
                MaxGap = atoi(optarg);
                break;

            case 'O':
                outAlignLen = atoi(optarg);
                break;

            case 'G':
                max_sub_rt_gap = atoi(optarg);
                break;

            case 'T':
                min_sub_rt_count = atoi(optarg);
                break;

            case 'r':
                PBS_minLen = atoi(optarg);
                break;

            case 'j':
                JoinThreshold = atof(optarg);
                break;

            case 'J':
                SplitThreshold = atof(optarg);
                break;

            case 'S':
                minOutScore = atof(optarg);
                break;

            case 'B':
                HigherSharpness = atof(optarg);
                break;

            case 'b':
                LowerSharpness = atof(optarg);
                break;

            case 'M':
                minMatchSim = atof(optarg);
                break;

            case 's':
                tRNA_file = optarg;
                break;

            case 'a':
                ps_dir = optarg;
                break;

            case 'f':
                fig_file = optarg;
                break;

            case 'C':
                checkCentriole = true;
                break;

            default:
                display_usage = true;

                //          case 'm': show_score_matrix(); break;
        }

    }

    if (optind >= argc || display_usage) {
        fprintf(stderr, "ltr_finder v%s\n", version);
        fprintf(stderr, "Usage  : [options] <INPUT_FASTA_FILE>\n");
        //  fprintf (stderr, "Options: -b NUM     bandwidth, default is %d\n",
        //     band_width);
        fprintf(stderr,
                "         -o NUM     gap open penalty, default is %d\n",
                gap_open);
        fprintf(stderr,
                "         -t NUM     gap extension penalty, default is %d\n",
                gap_ext);
        fprintf(stderr,
                "         -e NUM     gap end penalty, default is %d\n",
                gap_end);
        fprintf(stderr, "         -m NUM     match score, default is %d\n",
                score_match);
        fprintf(stderr, "         -u NUM     unmatch score, default is %d\n",
                score_mismatch);
        fprintf(stderr,
                "         -D NUM     Max distance between 5'&3'LTR, default is %d\n",
                Dmax);
        fprintf(stderr,
                "         -d NUM     Min distance between 5'&3'LTR, default is %d\n",
                Dmin);
        fprintf(stderr,
                "         -L NUM     Max length of 5'&3'LTR, default is %d\n",
                Lmax);
        fprintf(stderr,
                "         -l NUM     Min length of 5'&3'LTR, default is %d\n",
                Lmin);
        fprintf(stderr,
                "         -p NUM     min length of exact match pair, default is %d\n",
                Lex);
        fprintf(stderr,
                "         -g NUM     Max gap between joined pairs, default is %d\n",
                MaxGap);
        fprintf(stderr,
                "         -G NUM     Max gap between RT sub-domains, default is %d\n",
                max_sub_rt_gap);
        fprintf(stderr,
                "         -T NUM     Min sub-domains found in a RT domain, default is %d\n",
                min_sub_rt_count);
        fprintf(stderr, "         -j NUM     Threshold for join new sequence in existed alignment\n");
        fprintf(stderr, "                    new alignment similarity higher than this will be joined,\n");
        fprintf(stderr, "                    default is %0.2f\n",
                JoinThreshold);
        fprintf(stderr, "         -J NUM     Threshold for split existed alignment to two part\n");
        fprintf(stderr, "                    new alignment similarity lower than this will be split,\n");
        fprintf(stderr, "                    set this threshold lower than -j, means turn it off,\n");
        fprintf(stderr, "                    default is %0.2f\n",
                SplitThreshold);

        fprintf(stderr,
                "         -S NUM     output Score limit, default is %0.2f, [0,10]\n",
                minOutScore);
        fprintf(stderr,
                "         -M NUM     min LTR similarity threshold, default is %0.2f, [0,1]\n",
                minMatchSim);
        fprintf(stderr, "         -B NUM     Boundary alignment sharpness threshold, higher one.\n");
        fprintf(stderr, "                     one of the two edge's sharpness must higher than\n");
        fprintf(stderr, "                     this threshold, default is %0.3f, [0,1]\n",
                HigherSharpness);
        fprintf(stderr, "         -b NUM     Boundary alignment sharpness threshold, lower one.\n");
        fprintf(stderr, "                     both of the two edge's sharpness must higher than\n");
        fprintf(stderr, "                     this threshold, default is %0.3f, [0,1]\n",
                LowerSharpness);

        fprintf(stderr,
                "         -r NUM     PBS detecting threshold, min tRNA match length: %d, [1,18]\n",
                PBS_minLen);
        fprintf(stderr,
                "         -w NUM     output format: [0]-full, 1-summary, 2-table.\n");
        fprintf(stderr,
                "         -O NUM     output alignment length(only affect -w0), default is %d\n",
                outAlignLen);
        fprintf(stderr, "         -P STR     SeqIDs, will only calculate matched SeqID\n");
        fprintf(stderr, "                      POSIX style regular express is supported.\n");
        fprintf(stderr,
                "         -s filename      tRNA sequence file(FASTA format)\n");
        fprintf(stderr,
                "         -f filename      data file used to draw figure\n");
        fprintf(stderr,
                "         -a ps_scan_dir   Use ps_scan to predict protein domain\n");
        fprintf(stderr, "         -x         Output in html format\n");
        fprintf(stderr, "         -E         LTR must have edge signal\n");
        fprintf(stderr, "                    (at least two of PBS,PPT,TSR)\n");
        fprintf(stderr, "         -C         detect Centriole, delete highly repeat regions\n");
        fprintf(stderr,
                "         -F 01string      Filter to choose desired result,default is 0\n");
        fprintf(stderr, "                     10000000000 5'-LTR must have TG\n");
        fprintf(stderr, "                     01000000000 5'-LTR must have CA\n");
        fprintf(stderr, "                     00100000000 3'-LTR must have TG\n");
        fprintf(stderr, "                     00010000000 3'-LTR must have CA\n");
        fprintf(stderr, "                     00001000000 TSR must be found\n");
        fprintf(stderr, "                     00000100000 PBS must be found\n");
        fprintf(stderr, "                     00000010000 PPT must be found\n");
        fprintf(stderr, "                     00000001000 RT domain muse be found\n");
        fprintf(stderr, "                     00000000100 Integrase core must be found\n");
        fprintf(stderr, "                     00000000010 Integrase c-term must be found\n");
        fprintf(stderr, "                     00000000001 RNase H must be found\n");

        //  fprintf(stderr, "         -m         show score matrix\n");
        fprintf(stderr, "         -h         help\n");
        exit(1);
    }
    //    if(argc<2)
    //    {
    //        cerr<<"ltr_finder INPUT_FASTA_FILE"<<endl;
    //        exit(1);
    //    }

    //Program    : BGF
    //Version    : 2.1.2
    //Time       : Wed Nov 22 20:53:23 2006

    if (ishtml)
        printf("<html>\n<head>\n<title>LTR_FINDER Result</title>\n</head>\n<body>\n<pre>\n");

    printf("Program    : LTR_FINDER\n");

    printf("Version    : %s\n\n", version);

    const char *id_filter_stat;

    re_syntax_options = RE_SYNTAX_POSIX_EGREP |
            RE_BACKSLASH_ESCAPE_IN_LISTS | RE_DOT_NOT_NULL;

    id_filter_stat = re_compile_pattern(namePattern.c_str(),
            namePattern.length(), &id_filter);

    if (id_filter_stat != NULL) {
        printf("not a vaild POSIX regex after -P, code = %s\n", id_filter_stat);
    }

    if (tRNA_file != NULL) {
        timehere.markbeg();
        string tmp_out = "Load tRNA db [";
        tmp_out += tRNA_file;
        tmp_out += "] ";
        pbs.LoadSeq(tRNA_file);
        timehere.markend();
        timehere.outtime(tmp_out.c_str());
    }

    ps_scan.init(ps_dir);
    set_score_matrix(gap_open, gap_ext, gap_end, score_match, score_mismatch);

    FILE * inFASTA = fopen(argv[optind], "r");

    if (inFASTA == NULL) {
        cerr << "open " << argv[optind] << " error!" << endl;
        exit(1);
    }

    //init transDNA
    for (int i = 0; i < 128; ++i)
        transDNA[i] = 'N';
    transDNA['a'] = 'A';
    transDNA['A'] = 'A';
    transDNA['c'] = 'C';
    transDNA['C'] = 'C';
    transDNA['g'] = 'G';
    transDNA['G'] = 'G';
    transDNA['t'] = 'T';
    transDNA['T'] = 'T';
    transDNA['u'] = 'T';
    transDNA['U'] = 'T';

    seq_t sequence;
    sequence.s = NULL;
    sequence.m = 0;
    char name[1024];

    ofstream fig;

    if (fig_file != NULL)
        fig.open(fig_file);

    int total_sequence = 0;

    int total_img = 0;

    // read output of grf to sticks
    string grf_out = argv[optind + 1];
    map<string, vector <stick> > mst;
    ifstream in(grf_out.c_str());
    if (!in.good()) {
        cerr << "can't open " << grf_out << endl;
    }
    string line;
    while (getline(in, line)) {
        if (line[0] == '>') {
            vector<string> v;
            istringstream is(line);
            string s;
            while (getline(is, s, ':')) {
                v.push_back(s);
            }
            string chrom = v[0].substr(1);
            int start = atoi(v[1].c_str());
            int end = atoi(v[2].c_str());
            string cigar = v[3];
            
            int mismatch = 0, insertion = 0, deletion = 0, match = 0, count = 0;
            for (size_t i = 0; i < cigar.size(); i++) {
                char diff = cigar[i] - '0';
                if (diff >= 0 && diff <= 9) {
                    // is a digit
                    count = count * 10 + diff;
                } else {
                    // is a character
                    if (cigar[i] == 'm') {
                        match += count;
                    } else if (cigar[i] == 'M') {
                        mismatch += count;
                    } else if (cigar[i] == 'I') {
                        insertion += count;
                    } else if (cigar[i] == 'D') {
                        deletion += count;
                    }
                    count = 0;
                }
            }

            stick tmp;
            tmp.pos1 = start;
            tmp.end1 = start + match + mismatch + insertion - 1;
            tmp.pos2 = end - match - mismatch - deletion + 1;
            tmp.end2 = end;
            tmp.match_len = insertion > deletion ? match + mismatch + insertion
                    : match + mismatch + deletion;
            tmp.match_score = match;
            mst[chrom].push_back(tmp);
        }
    }
    in.close();
    
    while (-1 != read_fasta(inFASTA, &sequence, name, 0)) {
        total_sequence++;
        
        if (namePattern.length() > 0 && id_filter_stat == NULL &&
                re_search(&id_filter, name, strlen(name),
                0, strlen(name), 0) < 0) //name Pattern not found
            //namePattern.find(name)==string::npos)
            continue;

        for (int i = 0; i < sequence.l; ++i)
            sequence.s[i] = transDNA[sequence.s[i]];

        myString = (char*) sequence.s;

        stringLength = sequence.l;

        if (showPairNum)
            cout << "Sequence:" << name << " Len:" << stringLength << " ";

//        //cout<<"begin GetPairs\n";
//        vector < candidate > Pair;
//
//        GetPairs(myString, stringLength, Lex, Lmax, Dmin, Dmax, Pair);
//
//        if (showPairNum)
//            continue;
//
//        //output pairs for checking
//        for (int i = 0; i < Pair.size() && CHECK_PAIRS; ++i) {
//            cout << setw(10);
//            cout << Pair[i].pos1 << "-" << Pair[i].pos1 + Pair[i].len - 1 << " * "
//                    << Pair[i].pos2 << "-" << Pair[i].pos2 + Pair[i].len - 1 << endl;
//        }
//        //constrcut sticks
//        //cout<<used[2]<<endl;
//
//        vector < stick > sticks; 
//
//        JoinPairs(myString, stringLength, Pair, sticks);
        
        vector <stick> & sticks = mst[name];

        //cerr<<"pairs after join:"<<sticks.size()<<endl;
        int count_do = 0;

        ps_scan.reset();     

        for (int i = 0; i < sticks.size(); ++i) {
            sticks[i].score = 0;
            //cout<<"begin extend pairs\n";

            if (!ExtendPairs(myString, stringLength, sticks[i]))
                continue;

            count_do++;

            //cout<<"begin findsignal\n";
            FindSignal(myString, stringLength, pbs, ppt, sticks[i]);

            ps_scan.AddRegion(sticks[i].end1 + 1, sticks[i].pos2 - 1);
        }
        
        timehere.markbeg();
        ps_scan.Predict(myString, stringLength);
        timehere.markend();
        timehere.outtime("Predict protein Domains");

        for (int i = 0; i < sticks.size(); ++i) {
            if (sticks[i].score < 0)
                continue;
            
            //cout<<"begin finddomain\n";
            ps_scan.Find(sticks[i].end1 + 1, sticks[i].pos2 - 1, sticks[i].motif);

            //cout<<"count score\n"; 
            CountScore(sticks[i]);

            //cerr<<"len:"<<sticks[i].match_len<<" score:"<<sticks[i].match_score<<endl;
            sticks[i].match_score = sticks[i].match_score / sticks[i].match_len; //score/average_len

            if (sticks[i].match_score > 1) {
                //cerr << "similarity bigger than1:" << sticks[i].match_score << endl;
                sticks[i].match_score = 1; //why???
            }

        }
        //cerr<<"after minSharpness filter:"<<count_do<<endl;
        
        //may be this function is needn't
        //EraseOverlap(sticks);

        stable_sort(sticks.begin(), sticks.end());
        //cout<<"out result\n";

        int result_count = OutPutResult(name, myString, sequence.l, sticks, fig);

        if (CHECK_PAIRS)
            ps_scan.PrintNoUsed();

        cout << endl;

        if (result_count && ishtml && fig.good()) {
            total_img++;
            cout << "</pre>\n";
            cout << "<img src=\"" << total_img << ".png\">\n";
            cout << "<pre>\n";
        }

    }

    if (total_sequence == 0) {
        cout << "No sequence found, please input FASTA format sequence and try again\n";
    }

    fig.close();
    fclose(inFASTA);
    free(sequence.s);
    timeall.markend();
    timeall.outtime("Total consume");

    if (ishtml)
        cout << "</pre>\n</body>\n</html>\n";

    return (0);
}
