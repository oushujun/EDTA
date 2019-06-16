// Copyright (C) 2016 shij@miamioh.edu

#ifndef DETECTDR_H
#define DETECTDR_H

#include <map>
#include <iostream>
#include "parameter.h"
#include "DetectIR.h"

using std::string;
using std::vector;
using std::map;
using std::ofstream;

class DetectDR : public DetectIR {
public:
    DetectDR(parameter param);
    void getStem(string & s, int start2, string & cigar, 
            unsigned & tr1, unsigned & tr2);
    void getStem2(string & s, int start2, string & cigar, 
            unsigned & tr1, unsigned & tr2);
    void detectStr(int i, vector<pair<int, int> > & S_TR_SumScore,
            const string & seq, vector<repeat> & v) override;
    int findStartError(const string & s, int seed_len, string & left, 
            int & error, int start2);
    void alignment(string &s, int start2, int & sa, int & sb, int & mla, 
            int & mlb, int & m, int & i, string & new_left);
};

#endif /* DETECTDR_H */
