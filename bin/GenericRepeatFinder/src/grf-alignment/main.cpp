// Copyright (C) 2016 shij@miamioh.edu

#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>

using namespace std;

string decompress2(const string & s) {
    int count = 0;
    string r = "";
    for (auto i : s) {
        char diff = i - '0';
        if (diff >= 0 && diff <= 9) {
            // is a digit
            count = count * 10 + diff;
        } else {
            // is a character
            r += string(count, i);
            count = 0;
        }
    }
    return r;
}

string getAlignment(const string & s, const string & tir, const string & type) {
    string s1 = "", s2 = "", s3 = "";
    string tmp = decompress2(tir);
    if (type == "1") {        
        // pointers to ends of s
        int i = 0, j = s.size() - 1;
        for (auto k : tmp) {
            if (k == 'm') {
                s1 += s[i];
                s2 += s[j];
                i++;
                j--;
                s3 += '|';
            } else if (k == 'M') {
                s1 += s[i];
                s2 += s[j];
                i++;
                j--;
                s3 += ' ';                
            } else if (k == 'I') {
                s1 += s[i];
                s2 += '-';
                i++;
                s3 += ' ';
            } else if (k == 'D') {
                s2 += s[j];
                s1 += '-';
                j--;
                s3 += ' ';
            }
        }
    } else {
        // length of second tdr
        int tdr_len = 0;
        for (auto k : tmp) {
            if (k != 'I') {
                tdr_len++;
            }
        }
        int i = 0, j = s.size() - tdr_len;
        for (auto k : tmp) {
            if (k == 'm') {
                s1 += s[i];
                s2 += s[j];
                i++;
                j++;
                s3 += '|';
            } else if (k == 'M') {
                s1 += s[i];
                s2 += s[j];
                i++;
                j++;
                s3 += ' ';                
            } else if (k == 'I') {
                s1 += s[i];
                s2 += '-';
                i++;
                s3 += ' ';
            } else if (k == 'D') {
                s2 += s[j];
                s1 += '-';
                j++;
                s3 += ' ';
            }
        }        
    }
    return s1 + "\n" + s3 + "\n" + s2;
}

void alignment(const string & type, const string & input, 
        const string & output) {
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
    while (getline(in, line)) {
        if (line[0] == '>') {
            out << line << endl;
            vector<string> lines;
            string item;
            stringstream iss(line);
            while (getline(iss, item, ':')) {
                lines.push_back(item);
            }
            string tir = lines[3];
            getline(in, line);
            out << line << endl;            
            string r = getAlignment(line, tir, type);
            out << r << endl;
        } else {
            out << line << endl;
        }
    }
    in.close();
    out.close();
}

int main(int argc, char** argv) {
    string name = "grf-alignment";
    string help = 
            "This program shows the alignment of inverted/direct repeats.\n"
            "Usage: " + name + " <type> <input_fasta> <output>\n"
            "type <int>  input sequence type; 1: inverted repeats (or MITEs); 2: direct repeats.\n"
            ;
    if (argc == 1 || string(argv[1]) == "-h") {
        cout << help << endl;
        exit(0);
    } else if (argc < 4) {
        cerr << "Please input correct parameters." << endl;
        exit(1);
    }
    if (string(argv[1]) == "1" || string(argv[1]) == "2") {
        alignment(argv[1], argv[2], argv[3]);     
    } else {
        cerr << "<type> is not correct." << endl;
        exit(1);
    }       
    return 0;
}
