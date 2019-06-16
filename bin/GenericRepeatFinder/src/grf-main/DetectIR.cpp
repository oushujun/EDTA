// Copyright (C) 2016 shij@miamioh.edu

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <cctype>
#include <utility>
#include <numeric>
#include <omp.h>
#include <exception>
#include <regex>
#include <unordered_set>
#include "DetectIR.h"
#include "functions.h"

using namespace std;

DetectIR::DetectIR() {    
}

DetectIR::DetectIR(parameter param) : p(param) {
    // adjust parameters
    adjustParam();    
    string suffix;
    if (p.format == 0) {
        suffix = "fasta";
    } else {
        suffix = "id";
    }
    file1 = p.output + "/perfect." + suffix;
    file2 = p.output + "/perfect.spacer." + suffix;
    file3 = p.output + "/imperfect." + suffix;
}

void DetectIR::adjustParam() {
    if (!p.percent || (!p.max_indel && !p.max_mismatch)) {
        // perfect mapping
        p.seed_mismatch = 0;
        p.max_indel = 0;
        p.max_mismatch = 0;
        p.percent = 0;
    }
}

void DetectIR::run() {
    // read input data
    fastaread(chroms, genomeData, p.input);
    for (size_t i = 0; i < chroms.size(); i++) {
        chrom_pos[chroms[i]] = i;
    }
    out1.open(file1);
    out2.open(file2);
    out3.open(file3);
    if (!out1.good()) {
        cerr << "can't open file: " << file1 << endl;
        exit(1);
    }
    if (!out2.good()) {
        cerr << "can't open file: " << file2 << endl;
        exit(1);
    }
    if (!out3.good()) {
        cerr << "can't open file: " << file3 << endl;
        exit(1);
    }
    for (size_t i = 0; i < chroms.size(); i++) {
        // find candidates
        identifyCandidate(i);
        // remove redundant repeats and output unique repeats
        reduce(i);
    }
    out1.close();
    out2.close();
    out3.close();
    
    // filter results based on TR and spacer length
    if (p.max_stem != INT_MAX 
            || p.min_spacer_len != 0 || p.max_spacer_len != INT_MAX) {
        filterByLen(file1);
        filterByLen(file2);
        filterByLen(file3);
    }      
}

void DetectIR::identifyCandidate(int count) {
    time(&rawtime);
    cout << "calculate cumulative scores in chromosome " 
            << chroms[count] << ": " << ctime(&rawtime) << endl;
    string & seq = genomeData[count];
    // 1. Check sequence.
    size_t seqSize = seq.size();
    if (seqSize < (size_t) (2 * p.seed + p.min_space)) {
        cerr << "Input genome sequence size is short. Can't perform detection."
                << endl;
        return;
    }
    for (char & c : seq) {
        c = std::toupper(c);
    }
    // adjust max space according to input genome size
    if (seqSize - 2 * p.seed < (size_t) p.max_space) {
        max_space = seqSize - 2 * p.seed;
    } else {
        max_space = p.max_space;
    }
    // 2. Mapping.
    vector< pair<int, int> > cumScore(seqSize);
    mapVector(seq, cumScore);
    // 3. Calculate cumulative score.
    cumsum(cumScore);
    // 4. Possible terminal repeats' cumulative scores.
    vector< pair<int, int> > S_TR_SumScore(seqSize - p.seed + 1);
    S_TR_SumScore[0] = cumScore[p.seed - 1];
    for (size_t i = p.seed; i < seqSize; i++) {
        int first = cumScore[i].first - cumScore[i - p.seed].first;
        int second = cumScore[i].second - cumScore[i - p.seed].second;
        S_TR_SumScore[i - p.seed + 1] = std::make_pair(first, second);
    }
    cumScore.clear();
    cumScore.shrink_to_fit();
    // 5. Detecting candidates for different lengths
    candidates.resize(max_space - p.min_space + 1);
    time(&rawtime);
    cout << "find repeats in chromosome " << chroms[count] << ": " 
            << ctime(&rawtime) << endl;
#pragma omp parallel for   
    for (int i = p.min_space; i <= max_space; i++) {
        detectStr(i, S_TR_SumScore, seq, candidates[i - p.min_space]);
    }
}

void DetectIR::detectStr(int i, vector< pair<int, int> > & S_TR_SumScore,
        const string & seq, vector<repeat> & v) {   
    // numbers of paired terminal repeats with distance i
    size_t num = S_TR_SumScore.size() - i - p.seed;
    // length of repeat
    int l = i + 2 * p.seed;
    int max_absSum = p.seed_mismatch * 2;
    for (size_t j = 0; j < num; j++) {
        // absolute value of sum of scores of two terminal repeats
        int absSum = abs(S_TR_SumScore[j].first
                + S_TR_SumScore[j + p.seed + i].first) 
                + abs(S_TR_SumScore[j].second
                + S_TR_SumScore[j + p.seed + i].second);
        // j is start pos in seq
        if (absSum <= max_absSum) {
            size_t start = j;
            size_t end = start + l - 1;
            
            // verify base one by one
            // bases in boundary must match         
            if (isPair(seq[start], seq[end])) {
                int unpair = 0;
                for (int k = 1; k < p.seed; k++) {
                    if (!isPair(seq[start + k], seq[end - k])) {
                        unpair++;
                    }
                }
                if (unpair <= p.seed_mismatch) {
                    // extend seed region
                    string cigar = "";
                    unsigned tr1 = 0, tr2 = 0;
                    string candidateSeq = seq.substr(start, l);
                    
                    if (p.max_indel) {
                        // allow gap
                        getStem2(candidateSeq, cigar, tr1, tr2);
                    } else {
                        // not allow gap
                        getStem(candidateSeq, cigar, tr1, tr2);
                    }
                    if (cigar != "") {
                        repeat tmp;
                        tmp.start = j;
                        tmp.end = j + l - 1;
                        tmp.cigar = move(cigar);
                        tmp.tr1 = tr1;
                        tmp.tr2 = tr2;
                        v.push_back(move(tmp));
                    }
                }
            }
        }
    }
}

// get the stem without indel

void DetectIR::getStem(const string & s, string & cigar, 
            unsigned & tr1, unsigned & tr2) {
    int error = 0;
    // position of unpaired bases
    vector<int> pos;
    string left = "";
    int len = s.size();
    for (int i = 0, j = len - 1; i < j; i++, j--) {
        if (isPair(s[i], s[j])) {
            left += 'm';
        } else {
            error++;
            if (error > p.max_mismatch) {
                break;
            } else {
                left += 'M';
                pos.push_back(i);
            }
        }
    }
    // remove the Ms at the end
    int end;
    for (end = left.size() - 1; left[end] == 'M'; end--);
    len = end + 1;
    left = left.substr(0, len);
    if (p.percent == 100) {
        if (len >= p.min_stem) {
            cigar = compress(left);
            tr1 = tr2 = len;
        }
    } else {
        pos.push_back(len);
        for (int i = pos.size() - 1; i >= 0; i--) {
            if (pos[i] >= p.min_stem && i * 100 <= p.percent * pos[i]) {
                cigar = compress(left.substr(0, pos[i]));
                tr1 = tr2 = pos[i];
                return;
            }
        }
    }
}

// get stem with indel

void DetectIR::getStem2(const string & s, string & cigar, 
            unsigned & tr1, unsigned & tr2) {
    string left = "", new_left = "";
    // direct compare base to find the first unpaired base
    int error = 0;
    int len = s.size();
    // position of first unpaired base
    int start = findStartError(s, p.seed, left, error);
    // size of seq needing alignment
    int size = (start ? len - 2 * start : 0);
    if (size <= 2) {
        // no need for alignment
        // seq is fully paired or <=2 base in the middle
        int stem = (len - size) / 2;
        if (!(stem < p.min_stem || (stem == p.min_stem && left[stem - 1] == 'M')
                || p.percent * stem < error * 100)) {
            cigar = compress(left);
            tr1 = tr2 = stem;
        }
    } else {
        // start point of first half seq to align
        int sa = start;
        // start point of second half seq to align
        int sb = len - start - 1;
        int m = error;
        int i = 0;
        alignment(s, sa, sb, m, i, new_left);
        if (new_left == "") {
            // no alignment extension
            if (!(start < p.min_stem 
                    || (start == p.min_stem && left[start - 1] == 'M')
                    || p.percent * start < error * 100)) {
                cigar = compress(left);
                tr1 = tr2 = start;
            }
        } else {
            cigar = compress(left + new_left);
            tr1 = sa;
            tr2 = len - 1 - sb;
        }
    }
}

void DetectIR::alignment(const string & s, int & sa, int & sb, int & m, int & i,
        string & new_left) {
    vector< vector< pair<int, char> > > matrix;
    string a, b;
    int len = sb - sa + 1;
    if (len <= 2 * p.block) {
        a = s.substr(sa, len);
        b = reverseComplement(a);
        // score, pre
        matrix.resize(len + 1);
        // half matrix
        iniMatrix2(matrix, p.indel);
    } else {
        a = s.substr(sa, p.block);
        b = reverseComplement(s.substr(sb - p.block + 1, p.block));
        // score, pre
        matrix.resize(p.block + 1);
        // whole matrix
        iniMatrix(matrix, p.indel);
    }
    fillMatrix(matrix, a, b, p.match, p.mismatch, p.indel);
    string new_left2 = "";
    // find best alignment
    // len of aligned part of string a, b
    int la = 0, lb = 0;
    findPath(matrix, a, b, new_left2, sa, s.size() - 1 - sb, m, i, la, lb);
    if (new_left2 != "") {
        new_left += new_left2;
        sa += la;
        sb -= lb;
        if (len > 2 * p.block 
                && (la >= p.br * p.block || lb >= p.br * p.block)
                && sb > sa) {
            // best alignment end point in the outer part of matrix
            // from the end point to construct another block            
            alignment(s, sa, sb, m, i, new_left);
        }
    }
}

void DetectIR::findPath(vector<vector<pair<int, char > > > & matrix,
        const string & a, const string & b, string & left, int l1, int l2,
        int & mnum, int & inum, int & la, int & lb) {
    int max_score = 0;
    int x = 0, y = 0, mismatch = mnum, indel = inum;
    for (int i = 1; i < (int) matrix.size(); i++) {
        for (int j = 1; j < (int) matrix[i].size(); j++) {
            if (matrix[i][j].first > max_score
                    && matrix[i][j].second == 'd' && a[j - 1] == b[i - 1]
                    && l1 + j >= p.min_stem && l2 + i >= p.min_stem) {
                int m2 = 0, i2 = 0;
                findError(matrix, i, j, m2, i2, a, b);
                m2 += mnum;
                i2 += inum;
                // check error number, and percent
                if (m2 <= p.max_mismatch && i2 <= p.max_indel
                        && p.percent * (l1 + l2 + i + j) >= (2 * m2 + i2) * 100) {
                    max_score = matrix[i][j].first;
                    x = i;
                    y = j;
                    mismatch = m2;
                    indel = i2;
                }
            }
        }
    }
    la = y;
    lb = x;
    mnum = mismatch;
    inum = indel;
    while (x || y) {
        if (matrix[x][y].second == 'd') {
            if (a[y - 1] == b[x - 1]) {
                left = 'm' + left;
            } else {
                left = 'M' + left;
            }
            x--;
            y--;
        } else if (matrix[x][y].second == 'h') {
            left = 'I' + left;
            y--;
        } else if (matrix[x][y].second == 'v') {
            left = 'D' + left;
            x--;
        }
    }
}

int DetectIR::findStartError(const string & s, int seed_len, string & left, 
        int & error) {
    int start = 0;
    int len = s.size();
    // within seed region
    for (int i = 0, j = len - 1; i < seed_len; i++, j--) {
        if (isPair(s[i], s[j])) {
            left += 'm';
        } else {
            left += 'M';
            error++;
        }
    }
    // outside seed region
    for (int i = seed_len, j = len - 1 - seed_len; i <= j; i++, j--) {
        if (isPair(s[i], s[j])) {
            left += 'm';
        } else {
            start = i;
            break;
        }
    }
    return start;
}

void DetectIR::reduce(int i) {
    time(&rawtime);
    cout << "remove redundant sequences and print: " << ctime(&rawtime) << endl;
    // create index
    vector<vector<long long> > index(genomeData[i].size());
    bool flag;
    long long id  = 1;
    // from large size to small size
    for (int j = candidates.size() - 1; j >= 0; j--) {
        for (size_t k = 0; k < candidates[j].size(); k++) {
            unsigned s1 = candidates[j][k].start;
            unsigned e2 = candidates[j][k].end;
            unsigned l = candidates[j][k].tr1;
            unsigned r = candidates[j][k].tr2;
            unsigned e1 = s1 + l - 1;
            unsigned s2 = e2 - r + 1;
            flag = true;
            unordered_set<long long> tmp2(index[e1].begin(), index[e1].end());
            unordered_set<long long> tmp3(index[s2].begin(), index[s2].end());
            unordered_set<long long> tmp4(index[e2].begin(), index[e2].end());
            for (auto pos : index[s1]) {
                // seq is included in other long seq
                if (tmp2.count(pos) and tmp3.count(-pos) and tmp4.count(-pos)) {
                    flag = false;
                    break;
                }
            }
            if (flag) {
                // add to pos to index                
                for (size_t m = s1; m <= e1; m++) {
                    index[m].push_back(id);
                }
                for (size_t m = s2; m <= e2; m++) {
                    index[m].push_back(-id);
                }
                // print
               outputCandidate(candidates[j][k], i);
            }
            id++;
        }
        candidates[j].clear();
        candidates[j].shrink_to_fit();
    }
    candidates.clear();
    candidates.shrink_to_fit();
}

void DetectIR::outputCandidate(repeat & r, int i) {
    ofstream * out;
    unsigned size = r.end - r.start + 1;
    std::smatch s;
    if (regex_search(r.cigar, s, regex("^(\\d+)m$"))) {
        if (size == 2 * stoul(s[1].str())) {
            // perfect
            out = &out1;
        } else {
            // perfect with spacer
            out  = &out2;
        }
    } else {
        // imperfect
        out = &out3;
    }
    // filter
    if (p.r != -1) {
        unsigned spacer = size - r.tr1 - r.tr2;
        if (spacer > size * p.r) {
            return;
        }
    }
    *out << ">" << chroms[i] << ":"
            << r.start + 1 << ":"
            << r.end + 1 << ":"
            << r.cigar
            << endl;
    if (p.format == 0) {
        *out << genomeData[i].substr(r.start, size) << endl;
    }    
}

void DetectIR::filterByLen(const string & file) {
    string newFile = file + ".filterd";
    string cmd = p.program_path + "/grf-filter " 
        + to_string(p.min_stem) + " "
        + to_string(p.max_stem) + " "
        + to_string(p.min_spacer_len) + " "
        + to_string(p.max_spacer_len) + " "
        + file + " " + newFile;

    if (system(cmd.c_str())) {
        cerr << "failed to execute command " << cmd << endl;
        exit(1);
    }

    // remove unfiltered results
    remove(file.c_str());

    // rename the filtered results
    rename(newFile.c_str(), file.c_str());
}
