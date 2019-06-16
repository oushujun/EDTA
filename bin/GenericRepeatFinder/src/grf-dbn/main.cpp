// Copyright (C) 2016 shij@miamioh.edu

#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>

using namespace std;

void getDBN(string & s, const string & tir) {
    int count = 0;
    size_t m = 0, n = s.size() - 1;
    for (auto i : tir) {
        char diff = i - '0';
        if (diff >= 0 && diff <= 9) {
            // is a digit
            count = count * 10 + diff;
        } else {
            // is a character
            if (i == 'm') {
                for (int j = 0; j < count; j++) {
                    s[m++] = '(';
                    s[n--] = ')';
                }
            } else if (i == 'M') {
                m += count;
                n -= count;
            } else if (i == 'I') {
                m += count;
            } else if (i == 'D') {
                n -= count;
            }
            count = 0;
        }
    }
}

void outputDBN(const string & input, const string & output) {
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
    string tir = "";
    while (getline(in, line)) {
        if (line[0] == '>') {
            out << line << endl;
            vector<string> lines;
            string item;
            stringstream iss(line);
            while (getline(iss, item, ':')) {
                lines.push_back(item);
            }
            tir = lines[3];
            size_t l = stoul(lines[2]) - stoul(lines[1]) + 1;
            string dbn = string(l, '.');
            getDBN(dbn, tir);
            getline(in, line);
            out << line << endl;
            out << dbn << endl;
        } else {
            out << line << endl;
        }
    }
    in.close();
    out.close();
}

int main(int argc, char** argv) {
    string name = "grf-dbn";
    string help = 
            "This program converts inverted repeat/mite sequences in FASTA format to DBN format.\n"
            "Usage: " + name + " <input_fasta> <output_dbn>\n"
            ;
    if (argc == 1 || string(argv[1]) == "-h") {
        cout << help << endl;
        exit(0);
    } else if (argc < 3) {
        cerr << "Please input correct parameters." << endl;
        exit(1);
    }
    outputDBN(argv[1], argv[2]);    
    return 0;
}
