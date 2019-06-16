// Copyright (C) 2016 shij@miamioh.edu

#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

using namespace std;

void findCon(vector<string> & v1, vector<string> & v2, ofstream & out) {
    int base[] = {'A', 'T', 'C', 'G', 'N'};
    string con(v2[0].size(), ' ');
    for (size_t i = 0; i < con.size(); i++) {
        int freq[] = {0, 0, 0, 0, 0};
        for (size_t j = 0; j < v2.size(); j++) {
            char tmp = v2[j][i];
            if (tmp == 'A') {
                freq[0]++;
            } else if (tmp == 'T') {
                freq[1]++;
            } else if (tmp == 'C') {
                freq[2]++;
            } else if (tmp == 'G') {
                freq[3]++;
            } else {
                freq[4]++;
            }
        }
        // find base with largest frequency
        con[i] = base[max_element(freq, freq + 5) - freq];
    }
    out << con << endl;
    for (size_t i = 0; i < v1.size(); i++) {
        string aln(con.size(), ' ');
        out << v1[i] << endl;
        out << v2[i] << endl;
        string & s = v2[i];
        for (size_t j = 0; j < con.size(); j++) {
            if (s[j] == con[j] && con[j] != 'N') {
                aln[j] = 'm';
            } else {
                aln[j] = 'M';
            }
        }
        out << aln << endl;
    }
}

void alignment(const string & input, const string & output) {
    ifstream in(input);
    if (!in.good()) {
        cerr << "can't open file:" << input << endl;
        exit(1);
    }
    ofstream out(output);
    if (!out.good()) {
        cerr << "can't open file:" << output << endl;
        exit(1);
    }
    string line;
    // id, seq
    vector<string> v1, v2;
    getline(in, line);
    while (getline(in, line)) {
        if (line[0] == '-') {
            out << line << endl;            
            findCon(v1, v2, out);
            v1.clear();
            v2.clear();
        } else if (line[0] == '>') {
            // id
            v1.push_back(line);
            vector<string> lines;
            string item;
            stringstream iss(line);
            while (getline(iss, item, ':')) {
                lines.push_back(move(item));
            }
        } else {
            // fasta seq
            v2.push_back(line);
        }
    }
    findCon(v1, v2, out);
    in.close();
    out.close();
}

int main(int argc, char** argv) {
    string name = "grf-alignment2";
    string help = 
            "This program shows the alignments and consensus sequences of interspersed repeats.\n"
            "Usage: " + name + " <input> <output>\n"
            "The input file must be the output file of grf-intersperse (i.e., interspersed_repeat.out) and it must include repeat sequences.\n"
            ;
    if (argc == 1 || string(argv[1]) == "-h") {
        cout << help << endl;
        exit(0);
    } else if (argc < 3) {
        cerr << "Please input correct parameters." << endl;
        exit(1);
    }
    alignment(argv[1], argv[2]); 
    return 0;
}
