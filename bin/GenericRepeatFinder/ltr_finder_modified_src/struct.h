/*
* =====================================================================================
* 
*       Filename:  struct.h
* 
*    Description:  
* 
*        Version:  1.0
*        Created:  2006年10月20日 09时10分42秒 CST
*       Revision:  none
*       Compiler:  gcc
* 
*         Author:   (), 
*        Company:  
* 
* =====================================================================================
*/
#ifndef STRUCT_H_LTR_FINDER
#define STRUCT_H_LTR_FINDER
#include <vector>
#include <string>
#include "stdaln_interface.h"//interface to Liheng's Align Lib

using namespace std;

struct candidate
{
    int pos1;
    int pos2;
    int len;
};

struct MOTIF
{
    int begin;
    int end;
    vector<int> p_begin; //full-protein-begin
    vector<int> p_end; //full-protein-end
    int frame; //0,1,2
    int end_frame; //0,1,2
    bool forward; //y or n
    bool used; //
    string name;
    bool operator < (const MOTIF & a) const
    {
        return begin < a.begin;
    }
    bool operator == (const MOTIF & a) const
    {
        return (begin == a.begin &&
                end == a.end &&
                name == a.name);
    }
};

struct PBS_PPT
{
    bool PBS;
    bool _PBS;
    string PBS_name;
    string _PBS_name;
    int PBS_begin;
    int PBS_end;
    int PBS_match;
    int _PBS_begin;
    int _PBS_end;
    int _PBS_match;
    string pbs_str;
    string _pbs_str;

    bool PPT;
    bool _PPT;
    int PPT_begin;
    int PPT_count;
    int _PPT_begin;
    int _PPT_count;
    PBS_PPT(): PPT_count(0), _PPT_count(0), PBS_name(""), _PBS_name(""){
    };
};

struct stick
{
    vector < candidate > candi;
    int pos1;
    int pos2;
    int len1;
    int end1;
    int len2;
    int end2;
    float score;
    float match_score;
    float sharpness5;
    float sharpness3;
    int match_len;
    /*
        string PBS_name;
        int PBS_begin;
        int PBS_end;
        string PBS_str;
     
        int PPT_begin;
        int PPT_count;
    */
    int tg_pos1;
    int ca_pos1;
    int tg_pos2;
    int ca_pos2;
    int tsr_pos1;
    int tsr_pos2;
    int tsr_len;

    string LTR5;
    string LTR3;
    string strand_str;
    int strand;

    PBS_PPT pp;
    vector<MOTIF> motif;
    int status;
    int domain_score;

    /*
     * bool operator <  (const stick& a) const
     * {
     * if(offset == a.offset)
     * return pos1<a.pos1;
     * else
     * return offset<a.offset;
     * }
     */
    stick(): status(0), tsr_len(0), tg_pos1( -1), ca_pos1( -1), tg_pos2( -1), ca_pos2( -1), strand(0), score(0), match_score(0){
    };
    bool operator < (const stick & a) const
    {
        //        return score > a.score;

        if (pos1 == a.pos1)
            return end2 < a.end2;
        else
            return pos1 < a.pos1;
    }

};

struct TG_CA_TSR
{
    int tg_pos1;
    int ca_pos1;
    int tg_pos2;
    int ca_pos2;
    int tsr_pos1;
    int tsr_pos2;
    float signal_score;
    int score;
    int status;
    bool operator < (const TG_CA_TSR & a) const
    {
        if ( (a.tsr_pos1 != a.tg_pos1 && tsr_pos1 != tg_pos1)
                || (a.tsr_pos1 == a.tg_pos1 && tsr_pos1 == tg_pos1) ) //both has TSR, or both have not
        {
            return score > a.score; //return the pos near ori-boundary
        }
        else
            return (tg_pos1 -tsr_pos1) > (a.tg_pos1 - a.tsr_pos1); //then, who has TSR?
    }

};

struct TSR
{
    int lpos;
    int rpos;
    int len;
    short TG;
    short CA;
    bool operator < (const TSR & a) const
    {
        return ((TG * 80 + CA * 80 - (rpos - lpos)) >
                (a.TG * 80 + a.CA * 80 - (a.rpos - a.lpos)));
    }
};

struct FourPos
{
    int a;
    int b;
    int c;
    int d;
    bool operator < (const FourPos & fp) const
    {
        if (a != fp.a)
            return a < fp.a;
        else if (b != fp.b)
            return b < fp.b;
        else if (c != fp.c)
            return c < fp.c;
        else
            return d < fp.d;
    }
};

static bool pos1_sort (const candidate & a, const candidate & b)
{
    if (a.pos1 == b.pos1)
        return a.pos2 < b.pos2;
    else
        return a.pos1 < b.pos1;
};

static bool pos2_sort (const candidate & a, const candidate & b)
{
    if (a.pos2 == b.pos2)
        return a.pos1 < b.pos1;
    else
        return a.pos2 < b.pos2;
};

static bool stick_pos1_sort (const stick & a, const stick & b)
{
    if (a.pos1 == b.pos1)
        return a.end2 < b.end2;
    else
        return a.pos1 < b.pos1;
};
extern int Dmin ; // 100;
extern int Dmax ; // 15000;
extern int Lmin ; // 100;
extern int Lmax ; // 2000;
extern int Lex ; // 20;    //pre limit
extern int MaxGap ; // 200;   //between two pair
extern int outside_ext ; // 8;  //outside_extand bp
extern int inside_ext ; // 10;
extern int TSRmin ; // 4;
extern int TSRmax ; // 6;
extern int PBS_region;
extern int PBS_minLen;
extern int PPT_region;
extern int max_sub_rt_gap;
extern int min_sub_rt_count;
//extern float Sensitive;
extern int CHECK_PAIRS;
extern int wrought;
extern float minOutScore;
extern float minMatchSim;
extern float LowerSharpness;
extern float HigherSharpness;
extern float SplitThreshold;
extern float JoinThreshold;
extern int FModule;
extern int M5TG;
extern int M5CA;
extern int M3TG;
extern int M3CA;
extern int MTSR;
extern int MPBS;
extern int MPPT;
extern int MRT;
extern int MCORE;
extern int MCT;
extern int MRH;
extern int outAlignLen;
extern char *fig_file;
extern bool showPairNum;
extern bool checkCentriole;

extern char transDNA[128];
extern char transAA[128];
#endif

