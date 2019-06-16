/*
* =====================================================================================
* 
*       Filename:  CallPsScan.h
* 
*    Description:  system call ps_scan
* 
*        Version:  1.0
*        Created:  2006年11月05日 11时36分37秒 CST
*       Revision:  none
*       Compiler:  gcc
* 
*         Author:   (), 
*        Company:  
* 
* =====================================================================================
*/

#ifndef CALLPSSCAN_LTR_FINDER_H
#define CALLPSSCAN_LTR_FINDER_H
#include "struct.h"
#include "regex.h"
#include <vector>
#include <string>
#include <stdio.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <fstream>
#include <iostream>

using namespace std;

struct CRegion
{
    int first;
    int second;
    bool vaild;
    CRegion(int a = 0, int b = 0): vaild(true){
        first = a;
        second = b;
    };
    bool operator < (const CRegion& a) const
    {
        return first < a.first;
    }
};

struct DP_rec
{
    int pos; //current pos
    char p_stat[8]; //pre state
    int p_i[8]; //previous node
    int score[8]; //8 phase score, 7: no domain
    int c_i[8]; //point to MetaDomain
};

struct MetaDomain
{
    int begin;
    int end;
    char name;
    char phase;
    //DP_rec dp[8];
    bool operator < (const MetaDomain& a) const
    {
        return end < a.end;
    }

};
//struct MetaDomain
//{
//    int begin;
//    int begin2;
//    char domain[8];//1 2 3 4 5 6 7 meta domain
//    char frame[8];//their phase
//    char name;
//    char phase;
//    char name2;
//    char phase2;
//    char score;
//    bool operator < (const MetaDomain& a) const
//    {
//        return begin < a.begin;
//    }
//};

class CPSSCAN
{
    string path;
    string cmd;
    string aa_file;
    string motif_file;
    char c2i[128];
    bool ready;
    vector<MOTIF> motif;
    vector< CRegion > region;

    //for RT-domain
    vector<MetaDomain> md;
    vector< vector<int> > cluster_pos;
    int innerDist[7];
    int patLen[7];

    struct re_pattern_buffer pattern_buffer[7];
    char fastmap[7][1 << 8];

    int trans_seq(const char *nt, int len, char *aa, int is_trans = 1);
    void push_metadomain(const char* seq, int begin, int len, int frame);
    void join_metadomain();
    void push_motif(int len);
    bool has_metadomain(MetaDomain& domain, char name);

public:
    CPSSCAN(char* p){
        init(p);
    };
    CPSSCAN(){
        ready = false;
    };
    void init(char* p);
    ~CPSSCAN();
    bool Predict(const char* seq, int len);
    bool Find(int begin, int end, vector<MOTIF>& res);
    void PrintNoUsed(void);
    bool isReady();
    void AddRegion(int begin, int end);
    void reset();
};

#endif

