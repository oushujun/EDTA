// Copyright (C) 2016 shij@miamioh.edu

#ifndef DETECTMITE_H
#define DETECTMITE_H

#include <map>
#include "DetectIR.h"
#include "parameter.h"

using std::string;
using std::vector;
using std::ofstream;
using std::map;
using std::pair;

struct mite {
    unsigned int start;
    unsigned int end;
    unsigned int tr1;
    unsigned int tr2;
    string tsd;
    string tir;
};

class DetectMITE : public DetectIR {
public:
    DetectMITE(parameter param);
    void run() override;
    void identifyCandidate(int count) override;
    void detectStr(int i, vector<pair<int, int> > & S_TR_SumScore,
            const string & seq, vector<mite> & v);
    void outputCandidate();
private:
    vector<vector<vector<mite> > > candidates;
};

#endif /* DETECTMITE_H */
