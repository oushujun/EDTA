// Copyright (C) 2016 shij@miamioh.edu

#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>

using namespace std;

int check(const string & v) {
    int i;
    try {
        i = stoi(v);
    } catch (std::exception & e) {
        cerr << v << " is invalid." << endl;
        exit(1);
    }
    if (i < 0) {
        cerr << v << " is invalid." << endl;
        exit(1);
    }
    return i;
}

int getTirLen(const string & s, char c) {
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

void filter(int min_repeat, int max_repeat, int min_spacer, int max_spacer, 
        const string & input, const string & output) {
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
    bool flag = false;
    while (getline(in, line)) {
        if (line[0] == '-') {
            out << line << endl;
        } else if (line[0] == '>') {
            flag = false;
            vector<string> lines;
            string item;
            stringstream iss(line);
            while (getline(iss, item, ':')) {
                lines.push_back(move(item));
            }
            int len1 = getTirLen(lines[3], 'D');
            int len2 = getTirLen(lines[3], 'I');
            int len3 = stoi(lines[2]) - stoi(lines[1]) + 1 - len1 - len2;
            if (len1 >= min_repeat && len1 <= max_repeat
                    && len2 >= min_repeat && len2 <= max_repeat
                    && len3 >= min_spacer && len3 <= max_spacer) {
                out << line << endl;
                flag = true;
            }
        } else if (flag) {
            out << line << endl;
        }        
    }
    in.close();
    out.close();    
}

int main(int argc, char** argv) {
    string name = "grf-filter";
    string help = 
            "This program filters inverted repeats/MITEs/direct repeats "
            "according to given spacer and terminal repeat (TR) lengths.\n"
            "Usage: " + name + " <min_TR_len> <max_TR_len> <min_spacer_len>"
            " <max_spacer_len> <input_fasta> <output>\n"
            ;
    if (argc == 1 || string(argv[1]) == "-h") {
        cout << help << endl;
        exit(0);
    } else if (argc < 7) {
        cerr << "Please input correct parameters." << endl;
        exit(1);
    }
    int min_repeat = check(argv[1]);
    int max_repeat = check(argv[2]);
    int min_spacer = check(argv[3]);
    int max_spacer = check(argv[4]);
    filter(min_repeat, max_repeat, min_spacer, max_spacer, argv[5], argv[6]);
    return 0;
}
