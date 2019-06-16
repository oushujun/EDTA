// Copyright (C) 2016 shij@miamioh.edu

#ifndef DETECTINTERSPERSE_H
#define DETECTINTERSPERSE_H

#include <vector>
#include <map>
#include <unordered_map>
#include <iostream>
#include <tuple>
#include <utility>
#include "parameter.h"

using std::string;
using std::vector;
using std::map;
using std::unordered_map;
using std::pair;
using std::tuple;
using std::make_pair;

struct repeat {
    bool strand;
    bool flag = true;
    unsigned chrom;
    unsigned start;
    unsigned end;
};

struct group {
    // consensus seq
    string con;
    // position for each repeat copy
    vector<repeat> pos;    
};

class DetectIntersperse {
public:
    DetectIntersperse(const parameter & param);
    void run();
    void fastaread(vector<string> & id, vector<string> & seq,
            const string & file);
    void score(size_t i);
    void mapVector(const string & seq,
            vector<tuple<unsigned, unsigned, unsigned, unsigned> > & nseq);
    void cumsum(vector<tuple<unsigned, unsigned, unsigned, unsigned> > & v);
    void transform(const string & s,
        tuple<unsigned, unsigned, unsigned, unsigned> & t);
    void groupScore(unordered_map<string, unsigned> & index, size_t chrom);
    bool lowComplex(const string & s);
    bool filterSeq(const string & s);
    double gcContent(const string & s);
    double seqComplexity(const string & seq);
    void findCopy(vector<pair<size_t, size_t> > & v, vector<group> & r);
    void extendR(group & v);
    void extendL(group & v);
    string reverseComp(const string & s);
    char complement(char c);
    char findCon(char base[], int freq[], int len, double t);
    void filter(group & g);
    void output();
private:
    pair<size_t, size_t> mask = make_pair(10000, 0);
    char base[5] = {'A', 'T', 'C', 'G', 'N'};
    parameter p;
    vector<string> chroms;
    vector<string> genomeData;
    vector<tuple<unsigned, unsigned, unsigned, unsigned> > cumScore;
    // chrom, index
    vector<vector<pair<size_t, size_t> > > groupedScore;
    vector<vector<group> > result;
    vector<group> result2;
};

#endif /* DETECTINTERSPERSE_H */
