// Copyright (C) 2016 shij@miamioh.edu

#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <map>
#include <unordered_set>
#include <algorithm>

using namespace std;

void split(const string & s, vector<string> & lines, char sep){
    string item;
    stringstream iss(s);
    while (getline(iss, item, sep)) {
        lines.push_back(move(item));
    }
}

void getTirLen(const string & s, unsigned & l1, unsigned & l2) {
    int count = 0;
    for (auto i : s) {
        char diff = i - '0';
        if (diff >= 0 && diff <= 9) {
            // is a digit
            count = count * 10 + diff;
        } else {
            // is a character
            if (i != 'D') {
                l1 += count;
            }
            if (i != 'I') {
                l2 += count;
            }
            count = 0;
        }
    }
}

void outputOverlap(const string & input, const string & genome, 
        const string & output) {
    // genome length
    map<string, unsigned> lens;    
    ifstream in(genome);
    if (!in.good()) {
        cerr << "can't open file:" << genome << endl;
        exit(1);
    }
    string line, chrom = "";
    while (getline(in, line)) {
        if (line[0] == '>') {
            istringstream is(line.substr(1));
            is >> chrom;
            lens[chrom] = 0;
        } else {
            lens[chrom] += line.size();
        }
    }
    in.close();
    // tir
    unsigned s1 = 0, e1 = 0, s2 = 0, e2 = 0;
    string content;
    map<string, vector<tuple<unsigned, unsigned, unsigned, unsigned, string> > >
            tirs;
    in.open(input.c_str());
    while (getline(in, line)) {
        if (line[0] == '>') {
            content = line + "\n";
            vector<string> lines;
            split(line, lines, ':');            
            chrom = lines[0].substr(1);
            s1 = stoul(lines[1]) - 1;
            e2 = stoul(lines[2]) - 1;
            unsigned l1 = 0, l2 = 0;
            getTirLen(lines[3], l1, l2);
            e1 = s1 + l1 - 1;
            s2 = e2 - l2 + 1;
        } else if (line[0] != '-') {
            content += line;
            tirs[chrom].push_back(make_tuple(s1, e1, s2, e2, content));
        }
    }
    in.close();
    ofstream out(output);
    if (!out.good()) {
        cerr << "can't open file:" << output << endl;
        exit(1);
    }    
    // find nested tir
    for (auto it = lens.begin(); it != lens.end(); it++) {
        vector<vector<unsigned> > index(it->second);
        auto & v = tirs[it->first];
        vector<vector<string> > result;
        // sort tirs by length desc
        sort(v.begin(), v.end(), 
                [](tuple<unsigned, unsigned, unsigned, unsigned, string> & a,
                tuple<unsigned, unsigned, unsigned, unsigned, string> & b) {
                    unsigned l1 = get<3>(a) - get<0>(a);
                    unsigned l2 = get<3>(b) - get<0>(b);
                    return l1 > l2;
                });
        unsigned id = 0;
        for (size_t i = 0; i < v.size(); i++) {
            vector<string> lines;
            unsigned start = get<0>(v[i]);
            unsigned end = get<3>(v[i]);
            bool flag = false;
            for (auto & j : index[start]) {
                unordered_set<unsigned> tmp(index[end].begin(), 
                        index[end].end());
                if (tmp.count(j)) {
                    flag = true;
                    result[j].push_back(get<4>(v[i]));
                }
            }
            if (!flag) {
                // mark the spacer region
                unsigned start2 = get<1>(v[i]) + 1;
                unsigned end2 = get<2>(v[i]) - 1;
                for (size_t j = start2; j <= end2; j++) {
                    index[j].push_back(id);
                }
                result.resize(id + 1);
                result[id].push_back(get<4>(v[i]));
                id++;
            }
        }
        for (size_t i = 0; i < result.size(); i++) {
            if (result[i].size() > 1) {
                out << "----------------------------------------------" << endl;
                for (size_t j = 0; j < result[i].size(); j++) {
                    out << result[i][j] << endl;
                }                
            }
        }
    }
    out.close();
}

int main(int argc, char** argv) {
    string name = "grf-nest";
    string help =
            "This program finds nested TIRs/MITEs in a TIR/MTIE FASTA file.\n"
            "Usage: " + name + " <input_fasta> <genome_fasta> <output_fasta>\n"
            ;
    if (argc == 1 || string(argv[1]) == "-h") {
        cout << help << endl;
        exit(0);
    } else if (argc < 4) {
        cerr << "Please input correct parameters." << endl;
        exit(1);
    }
    outputOverlap(argv[1], argv[2], argv[3]);
    return 0;
}
