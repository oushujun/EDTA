// Copyright (C) 2016 shij@miamioh.edu

#ifndef FAMILYCLUSTER_H
#define FAMILYCLUSTER_H

#include <string>
#include <vector>
#include <map>

using std::string;
using std::map;
using std::vector;
using std::pair;

struct parameter {
    int flank = 50;
    int match = 1;
    int mismatch = 1;
    int indel = 2;
    int thread = 1;
    int copy = 3;
    string input = "";
    string genome = "";
    string output = "";
};

struct miteStr {
    unsigned int start;
    unsigned int end;
    string chrom;
    string tsd;
    string tir;
    string seq;
    string lf;
    string rf;
};

class FamilyCluster {
public:
    FamilyCluster(parameter param);
    void run();
private:
    parameter p;
    vector<string> chroms;
    map<string, int> chrom_pos;
    vector<string> genomeData;
    vector< vector<miteStr> > miteSet;
    vector<int> exemplar;
    void fastaread(vector<string> & id, vector<string> & seq, 
            const string & file);
    void processCandidate();
    void outputMiteSet();
    void outputExemplar();
    void flankSeqFilter(vector<miteStr> & v);
    void exemplarSelect(vector<miteStr> & v, int & e);
    void fillMatrix2(vector< vector< pair<int, char > > > & matrix,
        const string & a, const string & b, int match, int mismatch,
        int indel);
    int fillMatrix3(vector< vector<int> > & matrix,
        const string & a, const string & b, int match, int mismatch,
        int indel);
    int findPath(vector< vector< pair<int, char > > > & matrix,
        const string & a, const string & b);
    int swalign(
        int match,
        int mismatch,
        int indel,        
        const string & a,
        const string & b);
    int swalign2(
        int match,
        int mismatch,
        int indel,
        const string & a,
        const string & b);
};

#endif /* FAMILYCLUSTER_H */

