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
#include <regex>
#include <omp.h>
#include <exception>
#include "DetectMITE.h"
#include "functions.h"

using namespace std;

DetectMITE::DetectMITE(parameter param) {
    p = param;
    adjustParam();
}

void DetectMITE::run() {
    // read input data
    fastaread(chroms, genomeData, p.input);
    for (size_t i = 0; i < chroms.size(); i++) {
        chrom_pos[chroms[i]] = i;
    }
    candidates.resize(chroms.size());
    // find candidates
    for (size_t i = 0; i < chroms.size(); i++) {
        // find candidates
        identifyCandidate(i);   
    }
    // Output candidates
    outputCandidate();
    
    // filter results based on TR and spacer length
    if (p.max_stem != INT_MAX 
            || p.min_spacer_len != 0 || p.max_spacer_len != INT_MAX) {
        filterByLen(p.output + "/candidate.fasta");
    }   
}

void DetectMITE::identifyCandidate(int count) {
    time(&rawtime);
    cout << "calculate cumulative scores in chromosome " 
            << chroms[count] << ": " << ctime(&rawtime) << endl;
    string & seq = genomeData[count];
    // 1. Check sequence.
    size_t seqSize = seq.size();
    if (seqSize < (size_t) (2 * p.seed + p.max_space + 2 * p.max_tsd)) {
        cerr << "Input genome sequence size is short. Can't perform detection."
                << endl;
        return;
    }
    for (char & c : seq) {
        c = std::toupper(c);
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
    candidates[count].resize(p.max_space - p.min_space + 1);
    time(&rawtime);
    cout << "find candidates in chromosome " << chroms[count] << ": " 
            << ctime(&rawtime) << endl;
#pragma omp parallel for
    for (int i = p.min_space; i <= p.max_space; i++) {
        detectStr(i, S_TR_SumScore, seq, candidates[count][i - p.min_space]);
    }
}

void DetectMITE::detectStr(int i, vector< pair<int, int> > & S_TR_SumScore,
        const string & seq, vector<mite> & v) {
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
            size_t end = j + l -1;
            
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
                    string candidateSeq = seq.substr(start, l);
                    
                    // detect tsd
                    string tsd = detectTSD(seq, start, end, p.min_tsd,
                            p.max_tsd);
                    if (tsd != "") {
                        // extend seed region
                        string cigar = "";
                        unsigned tr1 = 0, tr2 = 0;
                        if (p.max_indel) {
                            // allow gap
                            getStem2(candidateSeq, cigar, tr1, tr2);
                        } else {
                            // not allow gap
                            getStem(candidateSeq, cigar, tr1, tr2);
                        }
                        // Filter out low complexity sequence.
                        if (cigar != "" 
                                && filterLowComplex(candidateSeq, cigar)) {
                            mite tmp;
                            tmp.start = j;
                            tmp.end = j + l - 1;
                            tmp.tir = move(cigar);
                            tmp.tsd = move(tsd);
                            tmp.tr1 = tr1;
                            tmp.tr2 = tr2;
                            v.push_back(move(tmp));
                        }
                    }
                }
            }
        }
    }
}

void DetectMITE::outputCandidate() {
    time(&rawtime);
    cout << "print: " << ctime(&rawtime) << endl;
    string file;
    if (p.format == 0) {
        file = p.output + "/candidate.fasta";
    } else {
        file = p.output + "/candidate.id";
    }
    ofstream out(file);
    if (!out.good()) {
        cerr << "can't open file: " << file << endl;
        exit(1);
    }
    // print results sorted by size
    for (int j = 0; j < p.max_space - p.min_space + 1; j++) {
        for (size_t i = 0; i < candidates.size(); i++) {
            for (size_t k = 0; k < candidates[i][j].size(); k++) {
                // filter
                if (p.r != -1) {
                    unsigned size = candidates[i][j][k].end 
                            - candidates[i][j][k].start + 1;
                    unsigned spacer = size - candidates[i][j][k].tr1 
                            - candidates[i][j][k].tr2;
                    if (spacer > size * p.r) {
                        continue;
                    }                 
                }                
                out << ">" << chroms[i] << ":"
                        << candidates[i][j][k].start + 1 << ":"
                        << candidates[i][j][k].end + 1 << ":"
                        << candidates[i][j][k].tir << ":"
                        << candidates[i][j][k].tsd
                        << endl;
                if (p.format == 0) {
                    out << genomeData[i].substr(candidates[i][j][k].start,
                            candidates[i][j][k].end - candidates[i][j][k].start
                            + 1) << endl;
                }
            }
        }
    }
    out.close();
}
