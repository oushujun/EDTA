// Copyright (C) 2016 shij@miamioh.edu

#include <iostream>
#include <vector>
#include <fstream>
#include <utility>
#include <regex>
#include <omp.h>
#include <exception>
#include <cmath>
#include <sstream>
#include "FamilyCluster.h"

using std::cout;
using std::cerr;
using std::endl;
using std::stoi;
using std::to_string;
using std::exception;
using std::string;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::move;
using std::max;
using std::pair;
using std::istringstream;

FamilyCluster::FamilyCluster(parameter param) : p(param) {
}

void FamilyCluster::run() {
    // read genome data
    fastaread(chroms, genomeData, p.genome);
    for (size_t i = 0; i < chroms.size(); i++) {
        chrom_pos[chroms[i]] = i;
    }
    // group candidate to sets
    processCandidate();
    // filter family members and select exemplar
    exemplar.resize(miteSet.size());
#pragma omp parallel for
    for (size_t i = 0; i < miteSet.size(); i++) {
        // filter out sequences flanked with similar sequences
        flankSeqFilter(miteSet[i]);      
        // select exemplar
        if (!miteSet[i].empty()) {
            exemplarSelect(miteSet[i], exemplar[i]);
        } else {
            exemplar[i] = -1;
        }
    }
    // output the mite set. Each mite group will be divided by dash line.
    outputMiteSet();
    // output the exemplars.
    outputExemplar();
}

// read fasta file to vector

void FamilyCluster::fastaread(vector<string> & id, vector<string> & seq, 
        const string & file) {
    ifstream fid(file);
    if (!fid.good()) {
        cerr << "can't open " << file << endl;
        exit(1);
    }
    string line;
    string header = "";
    while (getline(fid, line)) {
        if (line[0] == '>') {
            istringstream ss(line);
            ss >> header;
            header = header.substr(1);
            id.push_back(header);
            seq.push_back("");
        } else {
            seq.back() += line;
        }
    }
    fid.close();
}

void FamilyCluster::processCandidate() {
    ifstream in(p.input);
    if (!in.good()) {
        cerr << "can't open file: " << p.input << endl;
        exit(1);
    }
    vector<string> v;
    string line;
    while (getline(in, line)) {
        if (line[0] == '>') {
            // at least 3 copies
            if (v.size() >= unsigned(p.copy)) {
                vector<miteStr> set;
                for (auto & s : v) {
                    std::regex e(">(\\S+):(\\S+):(\\S+):(\\S+):(\\S+)\\.{3}");
                    std::smatch sm;
                    if (std::regex_search(s, sm, e)) {
                        miteStr tmp;
                        tmp.chrom = sm[1].str();
                        tmp.start = stoul(sm[2].str()) - 1;
                        tmp.end = stoul(sm[3].str()) - 1;
                        tmp.tir = sm[4].str();
                        tmp.tsd = sm[5].str();
                        const string & s = genomeData[chrom_pos[tmp.chrom]];
                        size_t len = s.size();
                        tmp.seq = s.substr(tmp.start, tmp.end - tmp.start + 1);
                        // extract flanking seq
                        if (tmp.start >= unsigned(p.flank)) {
                            tmp.lf = s.substr(tmp.start - p.flank, p.flank);
                        } else {
                            tmp.lf = s.substr(0, tmp.start);
                        }
                        if (tmp.end < len - unsigned(p.flank)) {
                            tmp.rf = s.substr(tmp.end + 1, p.flank);
                        } else {
                            tmp.rf = s.substr(tmp.end + 1, len - tmp.end - 1);
                        }
                        // add to candidate set
                        set.push_back(move(tmp));
                    } else {
                        cerr << "Input file is not in correct format." << endl;
                        exit(1);
                    }
                }
                miteSet.push_back(move(set));
            }
            v.clear();
        } else {
            v.push_back(line);
        }
    }
    in.close();
}

void FamilyCluster::outputMiteSet() {
    string file = p.output + "/miteSet.fasta";
    ofstream fid(file);
    if (!fid.good()) {
        cerr << "can't open file:" << file << endl;
        exit(1);
    }
    for (size_t i = 0; i < miteSet.size(); i++) {
        int copy = miteSet[i].size();
        if (copy) {
            fid << "-----------------------------------------------------------"
                    << endl;
            for (size_t j = 0; j < miteSet[i].size(); j++) {
                fid << ">" << miteSet[i][j].chrom << ":"
                        << miteSet[i][j].start + 1 << ":"
                        << miteSet[i][j].end + 1 << ":"
                        << miteSet[i][j].tir << ":"
                        << miteSet[i][j].tsd << endl;
                fid << miteSet[i][j].seq << endl;
            }
        }
    }
    fid.close();
}

void FamilyCluster::outputExemplar() {
    string file = p.output + "/mite.fasta";
    ofstream fid(file);
    if (!fid.good()) {
        cerr << "can't open file: " << file << endl;
        exit(1);
    }
    for (size_t i = 0; i < exemplar.size(); i++) {
        int index = exemplar[i];
        if (index != -1) {
            fid << ">" << miteSet[i][index].chrom << ":"
                    << miteSet[i][index].start + 1 << ":"
                    << miteSet[i][index].end + 1 << ":"
                    << miteSet[i][index].tir << ":"
                    << miteSet[i][index].tsd << ":"
                    << miteSet[i].size() << endl;
            fid << miteSet[i][index].seq << endl;
        }
    }
}

void FamilyCluster::flankSeqFilter(vector<miteStr> & v) {
    // status for each mite in a group; true: keep; false: delete
    vector<bool> v2(v.size(), true);
    // pairwise alignment
    for (size_t i = 0; i < v.size() - 1; i++) {
        if (v2[i]) {
            string seq1_lf = v[i].lf;
            string seq1_rf = v[i].rf;
            for (size_t j = i + 1; j < v.size(); j++) {
                if (v2[j]) {
                    string seq2_lf = v[j].lf;
                    string seq2_rf = v[j].rf;
                    int max_size1 = max(seq1_lf.size(), seq2_lf.size());
                    int max_size2 = max(seq1_rf.size(), seq2_rf.size());
                    // similarity = match number / the longest seq len
                    if ((double) swalign2(p.match, p.mismatch, p.indel,
                            seq1_lf, seq2_lf) / max_size1 > 0.5 
                            || (double) swalign2(p.match, p.mismatch, p.indel,
                            seq1_rf, seq2_rf) / max_size2 > 0.5)
                    {
                        // delete the second mite          
                        v2[j] = false;
                    }
                }
            }
        }
    }
    // number of different copy
    int diff = 0;
    for (auto i : v2) {
        if (i) {
            diff++;
        }
    }
    if (diff >= p.copy) {
        vector<miteStr> v3;
        for (size_t i = 0; i < v.size(); i++) {
            if (v2[i]) {
                v3.push_back(move(v[i]));
            }
        }
        v.clear();
        v = move(v3);
    } else {
        v.clear();
    }
}

void FamilyCluster::exemplarSelect(vector<miteStr> & v, int & e) {
    // pairwise alignment for every two mites in the group
    // calculate the similarities
    // find the one with max similarities
    size_t elementNum = v.size();
    // alignment score
    int similarities[elementNum][elementNum];
    for (size_t i = 0; i < elementNum - 1; i++) {
        similarities[i][i] = 0;
        for (size_t j = i + 1; j < elementNum; j++) {
            // alignment score
            similarities[i][j] = swalign(p.match, p.mismatch, p.indel,
                    v[i].seq, v[j].seq);
            similarities[j][i] = similarities[i][j];
        }
    }
    similarities[elementNum - 1][elementNum - 1] = 0;
    vector<int> score(elementNum);
    for (size_t i = 0; i < elementNum; i++) {
        int sum = 0;
        for (size_t j = 0; j < elementNum; j++) {
            sum += similarities[i][j];
        }
        score[i] = sum;
    }
    e = std::max_element(score.begin(), score.end())
            - score.begin();
}

// smith waterman algorithm
// return alignment score

int FamilyCluster::swalign(
        int match,
        int mismatch,
        int indel,
        const string & a,
        const string & b) {
    int score = 0;
    int row = b.size() + 1;
    int col = a.size() + 1;
    vector< vector<int> > matrix(row);
    // initialize score
    matrix[0].resize(col);
    for (int j = 0; j < col; j++) { 
        matrix[0][j] = 0;
    }
    for (int i = 1; i < row; i++) {
        matrix[i].resize(col);
        matrix[i][0] = 0;
    }
    // fill the score matrix
    score = fillMatrix3(matrix, a, b, match, mismatch, indel);
    return score;
}


// smith waterman algorithm
// return the number of matched bases

int FamilyCluster::swalign2(
        int match,
        int mismatch,
        int indel,
        const string & a,
        const string & b) {
    int matchNum = 0;
    int row = b.size() + 1;
    int col = a.size() + 1;
    vector< vector< pair<int, char> > > matrix(row);
    // initialize score
    matrix[0].resize(col);
    for (int j = 0; j < col; j++) { 
        matrix[0][j] = std::make_pair(0, 'n');
    }
    for (int i = 1; i < row; i++) {
        matrix[i].resize(col);
        matrix[i][0] = std::make_pair(0, 'n');
    }
    // fill the score matrix
    fillMatrix2(matrix, a, b, match, mismatch, indel);
    // find best alignment       
    matchNum = findPath(matrix, a, b);
    return matchNum;
}

void FamilyCluster::fillMatrix2(vector< vector< pair<int, char > > > & matrix,
        const string & a, const string & b, int match, int mismatch,
        int indel) {
    for (size_t i = 1; i < matrix.size(); i++) {
        for (size_t j = 1; j < matrix[i].size(); j++) {
            // d h v
            vector<int> v(3);
            vector<char> v2 = {'d', 'h', 'v'};
            if (a[j - 1] == b[i - 1]) {
                v[0] = matrix[i - 1][j - 1].first + match;
            } else {
                v[0] = matrix[i - 1][j - 1].first - mismatch;
            }
            v[1] = matrix[i][j - 1].first - indel;
            v[2] = matrix[i - 1][j].first - indel;
            auto k = std::max_element(v.begin(), v.end()) - v.begin();
            // local alignment
            if (v[k] < 0) {
                v[k] = 0;
            }
            matrix[i][j].first = v[k];
            matrix[i][j].second = v2[k];
        }
    }
}

int FamilyCluster::fillMatrix3(vector< vector<int> > & matrix,
        const string & a, const string & b, int match, int mismatch,
        int indel) {
    int max = 0;
    for (size_t i = 1; i < matrix.size(); i++) {
        for (size_t j = 1; j < matrix[i].size(); j++) {
            vector<int> v(3);
            if (a[j - 1] == b[i - 1]) {
                v[0] = matrix[i - 1][j - 1] + match;
            } else {
                v[0] = matrix[i - 1][j - 1] - mismatch;
            }
            v[1] = matrix[i][j - 1] - indel;
            v[2] = matrix[i - 1][j] - indel;
            auto k = std::max_element(v.begin(), v.end()) - v.begin();
            // local alignment
            if (v[k] < 0) {
                v[k] = 0;
            }
            matrix[i][j]= v[k];
            if (matrix[i][j] > max) {
                max = matrix[i][j];
            }
        }
    }
    return max;
}

int FamilyCluster::findPath(vector< vector< pair<int, char > > > & matrix,
        const string & a, const string & b) {
    int matchNum = 0;
    int max_score = 0;
    int x = 0, y = 0;
    for (size_t i = 1; i < matrix.size(); i++) {
        for (size_t j = 1; j < matrix[i].size(); j++) {
            if (matrix[i][j].first > max_score) {
                max_score = matrix[i][j].first;
                x = i;
                y = j;
            }
        }
    }
    if (max_score) {
        while ((x || y) && matrix[x][y].first != 0) {
            if (matrix[x][y].second == 'd') {
                if (a[y - 1] == b[x - 1]) {
                    matchNum++;
                }
                x--;
                y--;
            } else if (matrix[x][y].second == 'h') {
                y--;
            } else if (matrix[x][y].second == 'v') {
                x--;
            }
        }
    }
    return matchNum;
}
