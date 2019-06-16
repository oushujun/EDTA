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
#include "DetectDR.h"
#include "functions.h"

using namespace std;

DetectDR::DetectDR(parameter param) : DetectIR(param) {    
}

void DetectDR::detectStr(int i, vector< pair<int, int> > & S_TR_SumScore,
        const string & seq, vector<repeat> & v) {
    size_t seq_len = seq.size();
    // numbers of paired terminal repeats with distance i
    size_t num = S_TR_SumScore.size() - i - p.seed;
    // length of repeat
    int l = i + 2 * p.seed;
    int max_absSum = p.seed_mismatch * 2;
    for (size_t j = 0; j < num; j++) {
        // absolute value of difference of scores of two terminal repeats
        int absSum = abs(S_TR_SumScore[j].first
                - S_TR_SumScore[j + p.seed + i].first) 
                + abs(S_TR_SumScore[j].second
                - S_TR_SumScore[j + p.seed + i].second);
        // j is start pos in seq
        if (absSum <= max_absSum) {
            // seq len with longest extension region
            int max_ext = l + i;
            size_t start = j;
            int len = max_ext;
            
            if (seq_len < j + max_ext) {
                len = seq_len - j;
            }
            
            // verify base one by one
            // bases in boundary must match         
            if (match(seq[start], seq[start + i + p.seed])) {
                int unpair = 0;
                for (int k = 1; k < p.seed; k++) {
                    if (!match(seq[start + k], seq[start + k + i + p.seed])) {
                        unpair++;
                    }
                }
                if (unpair <= p.seed_mismatch) {
                    // extend seed region
                    string cigar = "";
                    unsigned tr1 = 0, tr2 = 0;
                    string candidateSeq = seq.substr(start, len);
                    
                    if (p.max_indel) {
                        // allow gap
                        getStem2(candidateSeq, i + p.seed, cigar, tr1, tr2);
                    } else {
                        // not allow gap
                        getStem(candidateSeq, i + p.seed, cigar, tr1, tr2);
                    }
                    if (cigar != "") {
                        repeat tmp;
                        tmp.start = j;
                        tmp.end = j + candidateSeq.size() - 1;
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

void DetectDR::getStem(string & s, int start2, string & cigar, 
            unsigned & tr1, unsigned & tr2) {
    int error = 0;
    // position of unpaired bases
    vector<int> pos;
    string left = "";
    int len = s.size();
    // start2 is the start of second tdr
    for (int i = 0; i < len - start2; i++) {
        if (match(s[i], s[i + start2])) {
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
    int len2 = end + 1;
    left = left.substr(0, len2);
    if (p.percent == 100) {
        if (len2 >= p.min_stem) {
            s = s.substr(0, start2 + len2);
            cigar = compress(left);
            tr1 = tr2 = len2;
        }
    } else {
        pos.push_back(len2);
        for (int i = pos.size() - 1; i >= 0; i--) {
            if (pos[i] >= p.min_stem && i * 100 <= p.percent * pos[i]) {
                s = s.substr(0, start2 + pos[i]);
                cigar = compress(left.substr(0, pos[i]));
                tr1 = tr2 = pos[i];
            }
        }
    }
}

void DetectDR::getStem2(string & s, int start2, string & cigar, 
            unsigned & tr1, unsigned & tr2) {
    string left = "", new_left = "";
    // direct compare base to find the first unpaired base
    int error = 0;
    int len = s.size();
    // position of first unpaired base in first half
    int start = findStartError(s, p.seed, left, error, start2);
    // number of base left
    int size = (start ? len - start2 - start : 0);
    if (size <= 1) {
        // no need for alignment
        // second half of the seq is fully matched, or <= 1 base left
        int stem = len - start2 - size;
        if (!(stem < p.min_stem || (stem == p.min_stem && left[stem - 1] == 'M')
                || p.percent * stem < error * 100)) {
            s = s.substr(0, stem + start2);
            cigar = compress(left);
            tr1 = tr2 = stem;
        }
    } else {
        // start points for alignment
        int sa = start;
        int sb = start + start2;
        // max len of first and second half seqs left to be aligned
        int mla = start2 - sa;
        int mlb = len - sb;
        // current mismatch number
        int m = error;
        // current indel number
        int i = 0;
        // length of second half seq aligned
        alignment(s, start2, sa, sb, mla, mlb, m, i, new_left);
        if (new_left == "") {
            // no alignment extension
            // check current alignment
            if (!(start < p.min_stem 
                    || (start == p.min_stem && left[start - 1] == 'M')
                    || p.percent * start < error * 100)) {
                s = s.substr(0, start2 + start);
                cigar = compress(left);
                tr1 = tr2 = start;
            }
        } else {
            left += new_left;
            s = s.substr(0, sb);
            cigar = compress(left);
            tr1 = sa;
            tr2 = sb - start2;
        }
    }
}

int DetectDR::findStartError(const string & s, int seed_len, string & left, 
        int & error, int start2) {
    int start = 0;
    int len = s.size();
    // within seed region
    for (int i = 0; i < seed_len; i++) {
        if (match(s[i], s[i + start2])) {
            left += 'm';
        } else {
            left += 'M';
            error++;
        }
    }
    // outside seed region
    for (int i = seed_len; i < len - start2; i++) {
        if (match(s[i], s[i + start2])) {
            left += 'm';
        } else {
            start = i;
            break;
        }
    }
    return start;
}

void DetectDR::alignment(string &s, int start2, int & sa, int & sb, int & mla, 
        int & mlb, int & m, int & i, string & new_left) {
    int size = min(min(mla, mlb), p.block);
    string a = s.substr(sa, size);
    string b = s.substr(sb, size);
    replace(b, '?');
    // score, pre
    vector< vector< pair<int, char> > > matrix(size + 1);
    iniMatrix(matrix, p.indel);
    fillMatrix(matrix, a, b, p.match, p.mismatch, p.indel);
    string new_left2 = "";
    // find best alignment
    // len of aligned part of string a, b
    int la = 0, lb = 0;
    findPath(matrix, a, b, new_left2, sa, sb - start2, m, i, la, lb);
    if (new_left2 != "") {
        sa += la;
        sb += lb;
        mla = start2 - sa;
        mlb = s.size() - sb;
        new_left += new_left2;
        if (mla > 0 && mlb > 0 
                && (la >= p.br * p.block || lb >= p.br * p.block)) {
            // best alignment end point in the outer part of matrix
            // from the end point to start another matrix
            alignment(s, start2, sa, sb, mla, mlb, m, i, new_left);
        }
    }
}
