// Copyright (C) 2016 shij@miamioh.edu

#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <regex>
#include <utility>
#include <exception>
#include <algorithm>
#include <climits>
#include "functions.h"

using namespace std;

// read fasta file to vector

void fastaread(vector<string> & id, vector<string> & seq, const string & file) {
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

void iniMatrix(vector< vector< pair<int, char > > > & matrix, int indel) {
    int size = matrix.size();
    matrix[0].resize(size);
    matrix[0][0] = std::make_pair(0, 'n');
    for (int j = 1; j < size; j++) {
        matrix[0][j] = std::make_pair(-j * indel, 'h');
    }
    for (int i = 1; i < size; i++) {
        matrix[i].resize(size);
        matrix[i][0] = std::make_pair(-i * indel, 'v');
    }
}

void iniMatrix2(vector< vector< pair<int, char > > > & matrix, int indel) {
    int size = matrix.size();
    matrix[0].resize(size);
    matrix[0][0] = std::make_pair(0, 'n');
    for (int j = 1; j < size; j++) {
        matrix[0][j] = std::make_pair(-j * indel, 'h');
    }
    for (int i = 1; i < size; i++) {
        int size2 = size - i;
        matrix[i].resize(size2);
        matrix[i][0] = std::make_pair(-i * indel, 'v');
    }
}

string reverseComplement(const string & s) {
    int size = s.size();
    string r(size, '?');
    for (int i = 0; i < size; i++) {
        switch (s[i]) {
            case 'A':
                r[size - 1 - i] = 'T';
                break;
            case 'C':
                r[size - 1 - i] = 'G';
                break;
            case 'G':
                r[size - 1 - i] = 'C';
                break;
            case 'T':
                r[size - 1 - i] = 'A';
                break;
            default:
                // not determined
                r[size - 1 - i] = '?';
                break;
        }
    }
    return r;
}

string reverse(const string & s) {
    int size = s.size();
    string r(size, ' ');
    for (int i = 0; i < size; i++) {
        r[size - 1 - i] = s[i];
    }
    return r;
}

void fillMatrix(vector< vector< pair<int, char > > > & matrix,
        const string & a, const string & b, int match, int mismatch,
        int indel) {
    vector<int> v(3);
    char v2[3] = {'d', 'h', 'v'};
    
    for (size_t i = 1; i < matrix.size(); i++) {
        for (size_t j = 1; j < matrix[i].size(); j++) {   
            if (a[j - 1] == b[i - 1]) {
                v[0] = matrix[i - 1][j - 1].first + match;
            } else {
                v[0] = matrix[i - 1][j - 1].first - mismatch;
            }
            v[1] = matrix[i][j - 1].first - indel;
            v[2] = matrix[i - 1][j].first - indel;
            auto k = std::max_element(v.begin(), v.end()) - v.begin();
            matrix[i][j].first = v[k];
            matrix[i][j].second = v2[k];
        }
    }
}

bool isPair(char a, char b) {
    if ((a == 'A' && b == 'T')
            || (a == 'C' && b == 'G')
            || (a == 'G' && b == 'C')
            || (a == 'T' && b == 'A')) {
        return true;
    } else {
        return false;
    }
}

string compress(const string & s) {
    char prev = s[0];
    int num = 0;
    string r = "";
    for (auto i : s) {
        if (i == prev) {
            num++;
        } else {
            r += to_string(num) + prev;
            prev = i;
            num = 1;
        }
    }
    r += to_string(num) + prev;
    return r;
}

string decompress(const string & s, char c) {
    int count = 0;
    string r = "";
    for (auto i : s) {
        char diff = i - '0';
        if (diff >= 0 && diff <= 9) {
            // is a digit
            count = count * 10 + diff;
        } else {
            // is a character
            if (i == 'm') {
                if (c == 'l') {
                    r += string(count, '(');
                } else {
                    r = string(count, ')') + r;
                }
            } else if (i == 'M') {
                if (c == 'l') {
                    r += string(count, '.');
                } else {
                    r = string(count, '.') + r;
                }
            } else if (i == 'I') {
                if (c == 'l') {
                    r += string(count, '.');
                }
            } else if (i == 'D') {
                if (c == 'r') {
                    r = string(count, '.') + r;
                }
            }
            count = 0;
        }
    }
    return r;
}

double seqComplexity(const string & seq) {
    unsigned int c = 1;
    string S(1, seq[0]);
    string Q = "", SQ = "";
    size_t seqLen = seq.size();
    for (size_t i = 1; i < seqLen; i++) {
        Q += seq[i];
        SQ = S + Q;
        string SQv = SQ.substr(0, SQ.size() - 1);
        if (SQv.find(Q) == string::npos) {
            S = SQ;
            Q = "";
            c++;
        }
    }
    return (c / (seqLen / (log(seqLen) / log(4))));
}

void cumsum(vector< pair<int, int> > & v) {
    pair<int, int> sum = {0, 0};
    for (size_t i = 0; i < v.size(); i++) {
        sum.first += v[i].first;
        sum.second += v[i].second;
        v[i] = sum;
    }
}

// use vector to represent sequence
// A = (1, 0); C = (0, 1); G = (0, -1); T = (-1, 0)

void mapVector(const string & seq, vector< pair<int, int> > & nseq) {
    for (size_t i = 0; i < seq.size(); i++) {
        switch (seq[i]) {
            case 'A':
                nseq[i] = std::make_pair(1, 0);
                break;
            case 'C':
                nseq[i] = std::make_pair(0, 1);
                break;
            case 'G':
                nseq[i] = std::make_pair(0, -1);
                break;
            case 'T':
                nseq[i] = std::make_pair(-1, 0);
                break;
            default:
                nseq[i] = std::make_pair(100, 0);
                break;
        }
    }
}

void readFile(vector<string> & v, const string & file) {
    ifstream in(file);
    if (!in.good()) {
        cerr << "can't open file: " << file << endl;
        exit(1);
    }
    string line;
    while (getline(in, line)) {
        if (line[0] == '>') {
            string s = "";
            v.push_back(s);
        } else {
            v.back() += line;
        }
    }
    in.close();
}

int checkInt(const string & s, const string & v) {
    int i;
    try {
        i = stoi(v);
    } catch (std::exception & e) {
        cerr << s << " " << v << " is invalid." << endl;
        exit(1);
    }
    if (i <= 0) {
        cerr << s << " " << v << " is invalid." << endl;
        exit(1);
    }
    return i;
}

int checkInt2(const string & s, const string & v) {
    int i;
    try {
        i = stoi(v);
    } catch (std::exception & e) {
        cerr << s << " " << v << " is invalid." << endl;
        exit(1);
    }
    if (i < 0) {
        cerr << s << " " << v << " is invalid." << endl;
        exit(1);
    }
    return i;
}

float checkFloat(const string & s, const string & v) {
    float i;
    try {
        i = stof(v);
    } catch (std::exception & e) {
        cerr << s << " " << v << " is invalid." << endl;
        exit(1);
    }
    return i;
}

int getStemLen(const string & s, char c) {
    int num = 0;
    int count = 0;
    for (auto i : s) {
        char diff = i - '0';
        if (diff >= 0 && diff <= 9) {
            // is a digit
            count = count * 10 + diff;
        } else {
            // is a character
            if (i != c) {
                num += count;
            }
            count = 0;
        }
    }
    return num;
}

string detectTSD(const string & seq, unsigned int start, unsigned int end,
        int min, int max) {
    string tsd = "";
    for (int i = max; i >= min; i--) {
        // out of boundary
        if (start < unsigned(i) || end + i > seq.size() - 1) {
            continue;
        }
        if (seq.substr(start - i, i) == seq.substr(end + 1, i)) {
            tsd = seq.substr(start - i, i);
            break;
        }
    }
    return tsd;
}

bool filterLowComplex(const string & seq, const string & tir) {
    // gc or at content < 0.2
    int len1 = getStemLen(tir, 'D');
    int len2 = getStemLen(tir, 'I');
    string l = seq.substr(0, len1);
    string r = seq.substr(seq.size() - len2, len2);
    float gc1 = gcContent(l);
    float gc2 = gcContent(r);
    if (gc1 < 0.2 || gc1 > 0.8 || gc2 < 0.2 || gc2 > 0.8) {
        return false;
    }
    // homopolymer
    if (regex_search(l, regex("(\\w)\\1{7,}|(\\w\\w)\\2{3,}"))
            || regex_search(r, regex("(\\w)\\1{7,}|(\\w\\w)\\2{3,}"))) {
        return false;
    }
    // lzc < 0.675
    if (seqComplexity(seq) < 0.675) {
        return false;
    }
    return true;
}

float gcContent(const string & s) {
    float r = 0;
    for (auto i : s) {
        if (i == 'G' || i == 'C') {
            r++;
        }
    }
    return r / s.size();
}

bool match(char a, char b) {
    if (a == b && (a == 'A' || a == 'T' || a == 'G' || a == 'C')) {
        return true;
    } else {
        return false;
    }
}

void replace(string & s, char a) {
    for (auto & i : s) {
        if (i != 'A' && i != 'T' && i != 'C' &&i != 'G') {
            i = a;
        }
    }
}

void findError(vector< vector< pair<int, char > > > & matrix, int i, int j,
        int & m2, int & i2, const string & a, const string & b) {
    int x = i, y = j;
    while (x || y) {
        if (matrix[x][y].second == 'd') {
            if (a[y - 1] != b[x - 1]) {
                m2++;
            }
            x--;
            y--;
        } else if (matrix[x][y].second == 'h') {
            i2++;
            y--;
        } else if (matrix[x][y].second == 'v') {
            i2++;
            x--;
        }
    }
}
