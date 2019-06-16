// Copyright (C) 2016 shij@miamioh.edu

#ifndef DETECTIR_H
#define DETECTIR_H

#include <ctime>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include "parameter.h"

using std::string;
using std::vector;
using std::map;
using std::ofstream;
using std::pair;

struct repeat {
    unsigned int start;
    unsigned int end;
    unsigned int tr1;
    unsigned int tr2;
    string cigar;
};

class DetectIR {
public:
    DetectIR();
    DetectIR(parameter param);
    void adjustParam();
    void virtual run();
    void getStem(const string & s, string & cigar, 
            unsigned & tr1, unsigned & tr2);
    void getStem2(const string & s, string & cigar, 
            unsigned & tr1, unsigned & tr2);
    void virtual identifyCandidate(int count);
    void virtual detectStr(int i, vector<pair<int, int> > & S_TR_SumScore,
            const string & seq, vector<repeat> & v);
    void alignment(const string & s, int & sa, int & sb, int & m, int & i, 
            string & new_left);
    int findStartError(const string & s, int seed_len, string & left, 
            int & error);
    void findPath(vector<vector<pair<int, char > > > & matrix,
            const string & a, const string & b, string & left, int l1, int l2,
            int & mnum, int & inum, int & la, int & lb);
    void reduce(int i);
    void outputCandidate(repeat & r, int i);
    void filterByLen(const string & file);
protected:
    int max_space;
    parameter p;
    time_t rawtime;
    ofstream out1, out2, out3;
    string file1, file2, file3;
    vector<string> chroms;
    map<string, int> chrom_pos;
    vector<string> genomeData;
    vector<vector<repeat> > candidates;
};

#endif /* DETECTIR_H */
