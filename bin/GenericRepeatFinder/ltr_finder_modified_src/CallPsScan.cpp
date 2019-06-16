/*
* =====================================================================================
*
*       Filename:  CallPsScan.cpp
*
*    Description:  system call ps_scan
*
*        Version:  1.0
*        Created:  2006年11月05日 15时48分34秒 CST
*       Revision:  none
*       Compiler:  gcc
*
*         Author:   (), 
*        Company:  
*
* =====================================================================================
*/

#include "CallPsScan.h"
#include <assert.h>
#include <algorithm>
#include <string.h>
using namespace std;
/* trans_seq: translate nt to aa
 * original version written by LI Heng
 * modified by XU Zhao
 */
int CPSSCAN::trans_seq(const char *nt, int len, char *aa, int is_trans)
{
    int i;
    char *p, c1, c2, c3;
    const static char* const aln_trans_table_eu_here =
        "KKNNRRSSTTTTIMIIEEDDGGGGAAAAVVVVQQHHRRRRPPPPLLLL**YY*WCCSSSSLLFFX";
    /* 01234567890123456789012345678901234567890123456789012345678901234 */

    //assert(nt); assert(aa);
    //if (len%3 != 0) return 1; /* fail */
    p = aa;

    for (i = 0; i < len - 2; i += 3)
    {
        c1 = nt[i];
        c2 = nt[i + 1];
        c3 = nt[i + 2];

        if (c1 > 3 || c2 > 3 || c3 > 3)
        {
            *p++ = 64; /* X */
        }
        else
        {
            *p++ = (char)((c1 << 4) | (c2 << 2) | c3);
        }

    }

    if (is_trans)
    {
        int L = len / 3;

        for (i = 0; i < L; ++i)
            aa[i] = aln_trans_table_eu_here[(int)aa[i]];

        aa[L] = '\0';
    }

    return 0;
}

void CPSSCAN::init(char* p)
{
    //init pattern

    if (p == NULL)
        ready = false; //will not exec ps_scan
    else
    {
        ready = true;
        path = p;
    }

    for (int i = 0;i < 128;++i)
        c2i[i] = 0;

    c2i['A'] = 0;

    c2i['G'] = 1; //

    c2i['C'] = 2; //

    c2i['T'] = 3;

    aa_file = "/tmp/ltr_aa.txt.";

    motif_file = "/tmp/ltr_mt.txt.";

    //aa_file="./ltr_aa.txt.";
    //motif_file="./ltr_mt.txt.";
    char buf[1024];

    sprintf(buf, "%d%d", getpid(), time(0));

    aa_file += buf;

    motif_file += buf;

    cmd = "PATH=";

    cmd += path;

    cmd += ":/usr/bin ";

    //cmd +=" ps_scan.pl -p PS50878 -p PS50879 -p PS50994 -p PS51027 -o pff -d ";
    cmd += " ps_scan.pl -p PS50879 -p PS50994 -p PS51027 -o pff -d ";

    //                   RNaseH      IN core    IN c-term
    cmd += path;

    cmd += "/prosite.dat ";

    cmd += aa_file;

    cmd += "> ";

    cmd += motif_file;

    //compile pattern
    re_syntax_options = RE_SYNTAX_POSIX_EGREP | RE_BACKSLASH_ESCAPE_IN_LISTS | RE_DOT_NOT_NULL;

    static const char* regex_pattern[] = {
                                             "[A-Z]{14}[NCQHSTWY][A-Z][AILMFPV]{2}[A-Z][AILMFPV][A-Z][KNS][A-Z]{4}"
                                             , "[AILMFPV][YCR][A-Z]?[AILMFPV][A-Z]{25}"
                                             , "[A-Z]{6}[AILMFPV][A-Z]{2}[AILMFPV]D[AILMFPV][A-Z]{2}[AG][FY][A-Z]{2}[AILMFPV][A-Z]{15}"
                                             , "[A-Z]{8}[AILMFPV][A-Z]{2}[GRE][A-Z]{3}[NCQHSTWY][PVSGACT][A-Z][AILMFPV]{2}[A-Z]{3}[AILMFPV][A-Z]{11}"
                                             , "[AILMFPV][A-Z]{2}[YF][A-Z]DD[AILMFPV]{3}[A-Z]{7}"
                                             , "[A-Z]{15}[AILMFPV][A-Z][AILMFPV][A-Z]{2}[RDEHK]K[A-Z]{2}[AILMFPV]",
                                             "[A-Z]{5}[AILMFPV][LT]G[A-Z]{2}[AILMFPV]"};

    for (int i = 0;i < 7;++i)
    {
        pattern_buffer[i].allocated = 0;
        pattern_buffer[i].buffer = 0;
        pattern_buffer[i].fastmap = fastmap[i];
        pattern_buffer[i].translate = 0;

        const char *id;
        id = re_compile_pattern(regex_pattern[i], strlen(regex_pattern[i]), &(pattern_buffer[i]));

        if (id != NULL)
        {
            printf(" error on compiling regex1. code = %s\n", id);
            exit(1);
        }

        re_compile_fastmap(&(pattern_buffer[i]));
    }

    innerDist[0] = 0;
    innerDist[1] = 7;
    innerDist[2] = 33;
    innerDist[3] = 25;
    innerDist[4] = 13;
    innerDist[5] = 6;
    innerDist[6] = 22;
    patLen[0] = 26;
    patLen[1] = 29;
    patLen[2] = 34;
    patLen[3] = 35;
    patLen[4] = 17;
    patLen[5] = 25;
    patLen[6] = 11;
}
CPSSCAN::~CPSSCAN()
{
    //    for(int i=0;i<7;++i)
    //        regfree(&(pattern_buffer[i]));
}
bool CPSSCAN::isReady()
{
    return ready;
}
bool CPSSCAN::Predict(const char* seq, int len)
{
    int minAAlen = 150;
    //combine regions
    sort(region.begin(), region.end());

    for (int i = 0;i < (int)region.size() - 1;++i)
    {
        if (region[i].second >= region[i + 1].first)
        {
            region[i + 1].first =
                region[i].first < region[i + 1].first ? region[i].first : region[i + 1].first;
            region[i + 1].second =
                region[i].second > region[i + 1].second ? region[i].second : region[i + 1].second;
            region[i].vaild = false;
        }

    }
    //for(int i=0;i<region.size();++i)
    //{
    //    cout<<" region:"<<region[i].first<<","<<region[i].second<<" vaild:"<<region[i].vaild<<endl;
    //}

    char* aaSeq[6];
    char* numSeq;

    numSeq = new char[(sizeof(char) * (len + 1))];

    for (int i = 0;i < len;++i)
        numSeq[i] = c2i[toupper(seq[i])];

    numSeq[len] = '\0';

    ofstream of;

    if (ready)
    {
        of.open(aa_file.c_str());

        if (!of.good())
            return false;
    }

    int aa_len = len / 3;

    for (int i = 0;i < 6;++i)
        aaSeq[i] = new char[(sizeof(char) * (aa_len + 1))];

    for (int i = 0;i < 3;++i)
    {
        trans_seq(numSeq + i, len - i, aaSeq[i]);
        //of<<">F"<<i<<"_0 "<<aa_len<<endl;
        //of<<aaSeq[i]<<endl;

        for (int j = 0;j < region.size();++j)
        {
            if (region[j].vaild)
            {
                int aa_b = (int)((region[j].first - i) / 3);

                if (aa_b < 0)
                    aa_b = 0;

                int aa_l = (int)((region[j].second - region[j].first + 1) / 3);

                //cout<<"find metadomain in ["<<i<<"] "<<aa_b<<"-"<<aa_b+aa_l<<endl;
                push_metadomain(aaSeq[i], aa_b, aa_l, i);

                if (!ready)
                    continue;

                int aa_e = aa_b + aa_l - 1;

                int out_b = aa_b;

                int current_len = 0;

                for (int k = 0;k < (aa_e - aa_b + 1);++k)
                {
                    if (aaSeq[i][aa_b + k] == '*' || k == (aa_e - aa_b))
                    {
                        if (current_len > minAAlen) ///??? min aa len, in which to find domain
                        {
                            of << ">F" << i << "_" << out_b << " " << current_len << " aa_b:" << aa_b << " aa_e:" << aa_e << "\n";

                            for (int m = out_b;m < out_b + current_len;++m)
                                of << aaSeq[i][m];

                            of << endl;
                        }

                        current_len = 0;
                        out_b = aa_b + k + 1;
                    }
                    else
                        current_len++;
                }

            }

        }
        //fprintf(tmp_file,">F%d\n%s\n",i,aaSeq);
    }

    join_metadomain();
    push_motif(len);
    cluster_pos.clear();
    md.clear(); //forward end

    reverse(numSeq, numSeq + len);

    for (int i = 0;i < len;++i)
        numSeq[i] = 3 - numSeq[i];

    for (int i = 0;i < 3;++i)
    {
        trans_seq(numSeq + i, len - i, aaSeq[i + 3]);
        //of<<">R"<<i<<"\n"<<aaSeq[i+3]<<"\n";

        for (int j = 0;j < region.size();++j)
        {
            if (region[j].vaild)
            {
                int aa_b = (int)((len - region[j].second - 1) / 3);

                if (aa_b < 0)
                    aa_b = 0;

                int aa_l = (int)((region[j].second - region[j].first + 1) / 3);

                push_metadomain(aaSeq[i + 3], aa_b, aa_l, i + 3);

                if (!ready)
                    continue;

                int aa_e = aa_b + aa_l - 1;

                int out_b = aa_b;

                int current_len = 0;

                for (int k = 0;k < (aa_e - aa_b + 1);++k)
                {
                    if (aaSeq[i + 3][aa_b + k] == '*' || k == aa_e - aa_b)
                    {
                        if (current_len > minAAlen) ///??? min aa len, in which to find domain
                        {
                            of << ">R" << i << "_" << out_b << " " << current_len << "\n";

                            for (int m = out_b;m < out_b + current_len;++m)
                                of << aaSeq[i + 3][m];

                            of << endl;
                        }

                        current_len = 0;
                        out_b = aa_b + k + 1;
                    }
                    else
                        current_len++;
                }

            }

        }

    }
    join_metadomain();
    push_motif(len);
    cluster_pos.clear();
    md.clear(); //backward end

    if (ready) //push motif predicted by ps_scan
    {
        of.close();

        if (system(cmd.c_str()))
        {
            cerr << "execute ps_scan error\n";
            exit(1);
        }

        //string cmd2="cat ";
        //cmd2+=motif_file;
        //system(cmd2.c_str());
        ifstream mf(motif_file.c_str());

        if (!mf.good())
        {
            cerr << "ps_scan output not found\n";
            return false;
        }

        int begin, end;
        string sname, tmp, domain;

        while (mf >> sname >> begin >> end >> domain >> tmp >> tmp >> tmp >> tmp >> tmp)
        {

            MOTIF tif;
            tif.name = domain;

            if (sname[0] == 'F')
            {
                int base = sname[1] - '0';
                int aa_b = atoi(sname.substr(3).c_str());
                //cout<<"readin aa_b="<<aa_b<<" "<<sname.substr(2).c_str()<<endl;
                begin += aa_b;
                end += aa_b;
                tif.begin = (begin - 1) * 3 + base;
                tif.end = (end - 1) * 3 + base + 2;
                tif.forward = true;
                tif.p_begin.push_back(tif.begin);
                tif.p_end.push_back(tif.end);
                tif.frame = tif.end_frame = base;

            }
            else if (sname[0] == 'R')
            {
                int base = sname[1] - '0';
                int aa_b = atoi(sname.substr(3).c_str());
                begin += aa_b;
                end += aa_b;
                tif.end = len - ((begin - 1) * 3 + base) - 1; //make sure, begin<end
                tif.begin = len - ((end - 1) * 3 + base + 2) - 1;
                tif.forward = false;
                tif.p_begin.push_back(tif.begin);
                tif.p_end.push_back(tif.end);
                tif.frame = tif.end_frame = base;

            }
            else
            {
                fprintf(stderr, "wrong motif_seq_name:%s\n", sname.c_str());
                return false;
            }

            tif.used = false;
            motif.push_back(tif);
            //cout<<"name:"<<tif.name<<" F:"<<tif.forward<<endl;
            //cout<<" begin:"<<tif.begin<<" end:"<<tif.end<<endl;
        }

        mf.close();
        unlink(motif_file.c_str());
        unlink(aa_file.c_str());
    }


    for (int i = 0;i < motif.size();++i) //find possible ORF
    {

        if (motif[i].forward)
        {
            //tif.begin=(begin-1)*3+base;
            //tif.end=(end-1)*3+base+2;
            int begin = (motif[i].begin - motif[i].frame) / 3 + 1;
            int end = (motif[i].end - 2 - motif[i].frame) / 3 + 1;
            int base = motif[i].frame;

            for (int j = begin - 1;j >= 0;--j)
            {
                //cout<<aaSeq[base][j];

                if (aaSeq[base][j] == 'M')
                {
                    //cout<<endl;
                    motif[i].p_begin.push_back(j*3 + base);
                    //cout<<"   j:"<<j<<endl;
                    //cout<<"   possible begin:"<<motif[i].p_begin[motif[i].p_begin.size()-1]<<endl;
                }
                else if (aaSeq[base][j] == '*')
                {
                    motif[i].p_begin.push_back((j + 1)*3 + base);
                    break;
                }

            }
            base = motif[i].end_frame;

            for (int j = end;j < aa_len;++j)
            {
                if (aaSeq[base][j] == '*')
                {
                    motif[i].p_end.push_back(j*3 + base + 2);
                    //cout<<"   possible end:"<<motif[i].p_end[motif[i].p_end.size()-1]<<endl;
                    break;
                }

            }

        }
        else
        {
            //tif.end=len-((begin-1)*3+base)-1;//make sure, begin<end
            //tif.begin=len-((end-1)*3+base+2)-1;
            int begin = ( -motif[i].end + len - 1 - motif[i].frame) / 3 + 1;
            int end = ( -motif[i].begin + len - 1 - 2 - motif[i].frame) / 3 + 1;
            int base = motif[i].frame;

            for (int j = begin - 1;j >= 0;--j)
            {
                if (aaSeq[base + 3][j] == 'M')
                {
                    motif[i].p_end.push_back(len - (j*3 + base) - 1);
                }
                else if (aaSeq[base + 3][j] == '*')
                {
                    motif[i].p_end.push_back(len - ((j + 1)*3 + base) - 1);
                    break;
                }

            }
            base = motif[i].end_frame;

            for (int j = end;j < aa_len;++j)
            {
                if (aaSeq[base + 3][j] == '*')
                {
                    motif[i].p_begin.push_back(len - (j*3 + base + 2) - 1);
                    break;
                }

            }

        }
        sort(motif[i].p_begin.begin(), motif[i].p_begin.end());
        sort(motif[i].p_end.begin(), motif[i].p_end.end());
    }

    for (int i = 0;i < 6;++i)
        delete[] aaSeq[i];

    delete[] numSeq;

    sort(motif.begin(), motif.end());

    return true;
}

bool CPSSCAN::Find(int begin, int end, vector<MOTIF>& res)
{
    //if(end<begin)
    //    return false;

    int i = 0;
    int j = motif.size();
    int m = j / 2;
    bool notyet = true;

    while (notyet)
    {
        //cout<<"i,j:"<<i<<" "<<j<<" m:"<<m<<endl;
        //cout<<"current:"<<motif[m].begin<<" need:"<<begin<<"-"<<end<<endl;

        if (i >= j)
            break;
        else if (motif[m].begin >= begin && motif[m].begin <= end)
            notyet = false;
        else if (motif[m].begin < begin)
        {
            i = m + 1;
            m = (int)((i + j) * 0.5);
            //cout<<"new :"<<i<<","<<m<<","<<j<<endl;
        }
        else if (motif[m].begin > end)
        {
            j = m - 1;
            m = (int)((i + j) * 0.5 + 0.5);
            //cout<<"new :"<<i<<","<<m<<","<<j<<endl;
        }

    }

    if (!notyet) //got it
    {

        for (i = m;i >= 0;--i)
        {
            if (motif[i].begin < begin)
                break;
            else if (motif[i].end <= end)
            {
                res.push_back(motif[i]);
                motif[i].used = true;
            }

        }

        for (i = m + 1;i < motif.size();++i)
        {
            if (motif[i].end > end)
                break;
            else if (motif[i].begin >= begin)
            {
                res.push_back(motif[i]);
                motif[i].used = true;
            }

        }

    }
    sort(res.begin(), res.end());

    return true;
}

void CPSSCAN::PrintNoUsed(void)
{
    int flag = 0;

    for (int i = 0;i < motif.size();++i)
    {
        if (!motif[i].used)
        {
            if (flag == 0)
            {
                flag = 1;
                cout << "*************DOMAIN NO USED*****************" << endl;
            }

            cout << "Domain: " << motif[i].name << endl;
            cout << "  from: " << motif[i].begin + 1
            << " to: " << motif[i].end + 1 << endl;
        }

    }
    motif.clear();
}
void CPSSCAN::AddRegion(int begin, int end)
{
    //if(!ready)
    //    return;

    bool flag = false;

    for (int i = 0;i < region.size();++i)
    {
        if (region[i].first <= end &&
                region[i].second >= begin)
        {
            region[i].first =
                region[i].first < begin ? region[i].first : begin;
            region[i].second =
                region[i].second > end ? region[i].second : end;
            flag = true;
            break;
        }

    }

    if (!flag)
    {
        CRegion cr(begin, end);
        region.push_back( cr );
    }
}
void CPSSCAN::reset()
{
    region.clear();
    motif.clear();
}
void CPSSCAN::push_metadomain(const char* seq, int begin, int len, int frame)
{
    MetaDomain tmp;
    vector<CRegion> rg;
    CRegion tmpRG;
    tmp.phase = frame;
    //tmp.score=1;
    //tmp.frame[0]=(char)frame;
    //struct re_registers regs;

    for (int i = 2;i < 6;++i)
    {
        tmp.name = i;
        //tmp.domain[0]=(char)i;
        int n = -1;

        while ( (n = re_search( &(pattern_buffer[i]), seq + begin, len, n + 1, len, 0)) >= 0 )
        {
            tmp.begin = n + begin;
            tmp.end = tmp.begin + patLen[i] - 1;
            md.push_back(tmp);

            tmpRG.first = n + begin - 270; //260 -> 270
            tmpRG.second = n + begin + 200; //188 -> 200
            rg.push_back(tmpRG);
        }

    }
    sort(rg.begin(), rg.end());

    for (int i = 0;i < (int)rg.size() - 1;++i)
    {
        if (rg[i].second >= rg[i + 1].first)
        {
            rg[i + 1].first =
                rg[i].first < rg[i + 1].first ? rg[i].first : rg[i + 1].first;
            rg[i + 1].second =
                rg[i].second > rg[i + 1].second ? rg[i].second : rg[i + 1].second;
            rg[i].vaild = false;
        }

    }

    for (int i = 0;i < rg.size();++i)
    {
        if (!rg[i].vaild)
            continue;

        if (rg[i].first < 0)
            rg[i].first = 0;

        if (rg[i].second >= begin + len)
            rg[i].second = begin + len;

        int rg_len = rg[i].second - rg[i].first + 1;

        int domain_name[] = {0, 1, 6};

        for (int j = 0;j < 3;++j)
        {
            tmp.name = domain_name[j];
            int n = -1;

            while ( (n = re_search( &(pattern_buffer[tmp.name]), seq + rg[i].first, rg_len, n + 1, rg_len, 0)) >= 0 )
            {
                tmp.begin = n + rg[i].first;
                tmp.end = tmp.begin + patLen[tmp.name] - 1;
                md.push_back(tmp);
            }

        }

    }

}
void CPSSCAN::join_metadomain()
{
    if (md.size() < 2) //needn't do DP
        return ;

    sort(md.begin(), md.end());

    //for(int i=0;i<md.size();++i)
    //    cout<<(int)md[i].name<<" ["<<(int)md[i].phase<<"] "<<md[i].begin<<"-"<<md[i].end<<endl;
    //cout<<"++++++++++++++++++++\n";
    //

    vector<DP_rec> dp;

    dp.reserve(md.size()); //a little larger

    DP_rec tdp;

    tdp.pos = 0;

    for (int i = 0;i < 8;++i) //first dp
    {
        tdp.p_stat[i] = -1;
        tdp.p_i[i] = -1;
        tdp.score[i] = 0;
        tdp.c_i[i] = -1;
    }

    dp.push_back(tdp);

    for (int i = 0;i < md.size();++i)
    {
        int score = 0;
        int pstat = -1;
        int pi = -1;
        //compute this state
        //neighbor

        for (int n = 1;n <= max_sub_rt_gap+1;++n) //max miss=2
        {
            int pre_stat = md[i].name - n;

            if (pre_stat < 0)
                pre_stat = 7; //no domain

            int min = 0;

            int max = md[i].begin; //default for pre_phase7

            if (pre_stat != 7)
            {
                for (int j = 1;j < n;++j)
                    min += patLen[pre_stat + j];

                max = min;

                for (int j = 1;j <= n;++j)
                    max += innerDist[pre_stat + j];
            }
            for (int j = (int)dp.size() - 1;j >= 0;--j)
            {
                if (md[i].begin - dp[j].pos >= min)
                {
                    if ( (pre_stat != 7 && dp[j].c_i[pre_stat] != -1 &&
                            md[i].begin - md[dp[j].c_i[pre_stat]].end <= max) //pre is not nondomain, and dist satisfy
                            || pre_stat == 7) //pre is nondomain, dist don't care
                    {

                        if (score < dp[j].score[pre_stat])
                        {
                            score = dp[j].score[pre_stat];
                            pstat = pre_stat;
                            pi = j;
                        }

                    }

                    break;
                }

            }

        }
        tdp = dp[dp.size() - 1]; //copy last

        if (tdp.pos == md[i].end) //has this node
        {

            if (score + 1 > tdp.score[md[i].name])
            {
                tdp.score[md[i].name] = score + 1;
                tdp.c_i[md[i].name] = i;
                tdp.p_stat[md[i].name] = pstat;
                tdp.p_i[md[i].name] = pi;
            }

        }
        else
        {
            tdp.score[md[i].name] = score + 1;
            tdp.c_i[md[i].name] = i;
            tdp.p_stat[md[i].name] = pstat;
            tdp.p_i[md[i].name] = pi;
            tdp.pos = md[i].end;
            dp.push_back(tdp);
        }
        //find max pre stat for phase7
        score = -1;

        pstat = -1;

        pi = dp.size() - 2;

        for (int j = 0;j < 8;++j)
        {
            if (score < dp[pi].score[j])
            {
                pstat = j;
                score = dp[pi].score[j];
            }

        }

        if (score >= dp[dp.size() - 1].score[7])
        {
            dp[dp.size() - 1].score[7] = score;
            dp[dp.size() - 1].p_stat[7] = pstat;
            dp[dp.size() - 1].p_i[7] = pi;
        }
        //        cout<<"md["<<i<<"] name:"<<(int)md[i].name<<" "<<md[i].begin<<"-"<<md[i].end<<endl;
        //        cout<<"       curr index:"<<dp.size()-1<<"  pos:"<<dp[dp.size()-1].pos<<endl;
        //        for(int j=0;j<8;++j)
        //        {
        //            cout<<"    stat:"<<j<<"  pre_stat:"<<(int)dp[dp.size()-1].p_stat[j]
        //                <<" p_i:"<<dp[dp.size()-1].p_i[j]
        //                <<" score:"<<dp[dp.size()-1].score[j]
        //                <<" c_i:"<<dp[dp.size()-1].c_i[j]<<endl;
        //        }

    }
    //trace back
    int pi = dp.size() - 1;

    int pstat = -1;

    int score = -1;

    vector<int> pos;

    for (int i = 0;i < 8;++i)
    {
        if (score < dp[pi].score[i])
        {
            score = dp[pi].score[i];
            pstat = i;
        }

    }

    while (pstat != -1)
    {
        if ( dp[pi].c_i[pstat] != -1 ) //has node
        {
            pos.push_back(dp[pi].c_i[pstat]);

            if (dp[pi].p_stat[pstat] == 7) //first ele
            {
                cluster_pos.push_back(pos);
                pos.clear();
            }

        }
        int old_stat = pstat;
        pstat = dp[pi].p_stat[old_stat];
        pi = dp[pi].p_i[old_stat];
    }
    //    for(int i=(int)cluster_pos.size()-1;i>=0;--i)
    //    {
    //        cout<<"------\n";
    //        for(int j=(int)(cluster_pos[i].size())-1;j>=0;--j)
    //        {
    //            int p=cluster_pos[i][j];
    //            cout<<(int)md[p].name<<"\t "<<md[p].begin<<"-"<<md[p].end<<" phase:"<<(int)md[p].phase<<endl;
    //        }
    //    }
    //    cout<<"++++++++++++\n";

    //clear it after call join_metadomain
    //    cluster_pos.clear();

}
//void CPSSCAN::join_metadomain()
//{
//    sort(md.begin(),md.end());
//    for(int i=0;i<(int)md.size()-1;++i)
//    {
//        for(int j=i+1;j<md.size();++j)
//        {
//            int d= md[j].begin-md[i].begin2-patLen[md[i].name2];
//            if(d > 140)
//                break;
//            //cout <<"judge: "<<(int)md[i].name2<<" "<<(int)md[j].name<<" pos:"<<md[i].begin2<<" to"<<md[j].begin<<endl;
//            if(md[j].name-md[i].name2 == 1)//neighbor metadomain
//            {
//                //cout<<" dist:"<<d<<endl;
//                if(d < 0 || d > innerDist[md[j].name-1])
//                    continue;
//                //if(d > 34)//innerDist[md[j].name-1])
//                //    break;
//                md[i].begin2=md[j].begin;
//                md[i].name2=md[j].name;
//                md[i].domain[md[i].score]=md[j].name;
//                md[i].frame[md[i].score]=md[j].phase;
//                md[i].score ++;
//            }
//            else if(md[j].name-md[i].name2 == 2)//next to neighbor metadomain
//            {
//                //cout<<" dist:"<<d<<endl;
//                if(d < patLen[md[j].name-1]
//                        || d > innerDist[md[j].name-1]+innerDist[md[j].name-2]+patLen[md[j].name-1])
//                    continue;
//                //if(d > 102)
//                //    break;
//                md[i].begin2=md[j].begin;
//                md[i].name2=md[j].name;
//                md[i].domain[md[i].score]=md[j].name;
//                md[i].frame[md[i].score]=md[j].phase;
//                md[i].score ++;
//            }
//            else if(md[j].name-md[i].name2 == 3)//6 to 3, 7 to 4, or so
//            {
//                //cout<<" dist:"<<d<<endl;
//                if(d < patLen[md[j].name-1] + patLen[md[j].name-2]
//                        || d > (innerDist[md[j].name-1]+innerDist[md[j].name-2]+
//                               innerDist[md[j].name-3]+patLen[md[j].name-1]+patLen[md[j].name-2]) )
//                    continue;
//                md[i].begin2=md[j].begin;
//                md[i].name2=md[j].name;
//                md[i].domain[md[i].score]=md[j].name;
//                md[i].frame[md[i].score]=md[j].phase;
//                md[i].score ++;
//            }
//
//        }
//    }
//}
//bool CPSSCAN::has_metadomain(MetaDomain& domain, char name)
//{
//    for(int i=0;i<domain.score;++i)
//    {
//        if(domain.domain[i] == name)
//            return true;
//    }
//    return false;
//}
void CPSSCAN::push_motif(int len)
{
    for (int i = 0;i < cluster_pos.size();++i) //push motif
    {
        //cout<<i<<" score:"<<md[i].score<<" "<<md[i].name<<" - "<<md[i].name2<<endl;
        //cout<<" has 3:"<<has_metadomain(md[i],2)<<endl;
        //cout<<" has 5:"<<has_metadomain(md[i],4)<<endl;
        //cout<<" has 7:"<<has_metadomain(md[i],6)<<endl;

        if (cluster_pos[i].size() >= min_sub_rt_count )
            //&& has_metadomain(md[i],2)
            //&& has_metadomain(md[i],4)
            //&& has_metadomain(md[i],6))
        {
            //cout<<"FIND!\n";
            MOTIF tif;
            int last = cluster_pos[i][0];
            int first = cluster_pos[i][cluster_pos[i].size() - 1];

            if (md[first].phase < 3) //forward
            {
                int base = md[first].phase;
                tif.begin = md[first].begin * 3 + base;
                tif.end = md[last].end * 3 + base + 2;
                tif.forward = true;
                tif.p_begin.push_back(tif.begin);
                tif.p_end.push_back(tif.end);
                tif.frame = base;
                tif.end_frame = md[last].phase;

                for (int j = (int)cluster_pos[i].size() - 1;j >= 0;--j)
                {
                    tif.name += '0' + md[cluster_pos[i][j]].name + 1;
                    tif.name += '[';
                    tif.name += '0' + md[cluster_pos[i][j]].phase;
                    tif.name += ']';
                }

            }
            else
            {
                int base = md[first].phase - 3;
                tif.end = len - (md[first].begin * 3 + base) - 1; //make sure, begin<end
                tif.begin = len - (md[last].end * 3 + base + 2) - 1;
                tif.forward = false;
                tif.p_begin.push_back(tif.begin);
                tif.p_end.push_back(tif.end);
                tif.frame = base;
                tif.end_frame = md[last].phase - 3;

                for (int j = 0;j < cluster_pos[i].size();++j)
                {
                    tif.name += '0' + md[cluster_pos[i][j]].name + 1;
                    tif.name += '[';
                    tif.name += '0' + md[cluster_pos[i][j]].phase;
                    tif.name += ']';
                }


            }

            if (!CHECK_PAIRS)
                tif.name = "RT";
            else
                tif.name = "RT(" + tif.name + ")";

            tif.used = false;

            motif.push_back(tif);
        }

    }
}

