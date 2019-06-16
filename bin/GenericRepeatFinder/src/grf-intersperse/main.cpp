// Copyright (C) 2016 shij@miamioh.edu

#include <cstdlib>
#include <iostream>
#include <string>
#include <exception>
#include <omp.h>
#include "DetectIntersperse.h"

using namespace std;

unsigned checkInt(const string & s, const string & v) {
    unsigned i;
    try {
        i = stoul(v);
    } catch (exception & e) {
        cerr << s << " is invalid." << endl;
        exit(1);
    }
    return i;
}

void checkParameter(int argc, char** argv, parameter & p) {
    string name = "grf-intersperse";
    string help =
            "This program performs genome-wide interspersed repeat detection "
            "for input genome sequences.\n"
            "Usage: " + name + " [options]\n"
            "[mandatory parameters]\n"
            "-i <string>  Input genome sequence file in FASTA format.\n"
            "-o <string>  Output directory.\n"
            "\n"
            "[optional parameters]\n"
            "-t <int>  Number of threads used in this program; default = 1.\n"
            "-f <int>  Format for outputs; 0: consensus sequence with position "
            "and sequence of each repeat copy; 1: consensus sequence with "
            "position of each repeat copy only; default = 0.\n"
            "-c <int>  Minimum copy number of seeds in the genome; "
            "default = 3.\n"
            "-s <int>  Length of the seed region; default = 20; must >= 10.\n"
            "-n <int>  Maximum number of undetermined bases in the "
            "extension of the seed region (either direction); default = 1.\n"
            "-p <int>  Minimum identity percentage for a determined base in "
            "the consensus sequence; default = 80, which means that any "
            "determined base in the consensus sequence must be identical to at "
            "least 80% of the bases in the same position from all repeat "
            "copies.\n"
            "-m <int>  Maximum number of mismatches in repeat copies compared "
            "with the determined bases of the consensus sequence; repeat copies"
            " with mismatches > this value will be removed from the group; "
            "default = 2.\n"
            ;
    if (argc == 1 || string(argv[1]) == "-h") {
        cout << help << endl;
        exit(0);
    }
    for (int i = 1; i < argc; i += 2) {
        string s(argv[i]);
        if (i + 1 >= argc) {
            cerr << "Total parameter number is not correct." << endl;
            exit(1);
        }
        string v(argv[i + 1]);
        if (s == "-i") {
            p.input = v;
        } else if (s == "-o") {
            p.output = v;         
        } else if (s == "-s") {
            p.seed = checkInt(s, v);
            if (p.seed < 10) {
                cerr << s << " is invalid." << endl;
                exit(1);
            }
        } else if (s == "-t") {
            p.thread = checkInt(s, v);
            if (p.thread < 1) {
                cerr << s << " is invalid." << endl;
                exit(1);
            }
        } else if (s == "-f") {
            p.format = checkInt(s, v);
            if (p.format != 0 && p.format != 1) {
                cerr << s << " is invalid." << endl;
                exit(1);                
            }
        } else if (s == "-c") {
            p.copy = checkInt(s, v);
            if (p.copy < 1) {
                cerr << s << " is invalid." << endl;
                exit(1);
            }
        } else if (s == "-n") {    
            p.max_n = checkInt(s, v);
        } else if (s == "-p") {    
            p.min_con = checkInt(s, v);
            if (p.min_con < 1 || p.min_con > 100) {
                cerr << s << " is invalid." << endl;
                exit(1);
            }
            p.min_con /= 100;
        } else if (s == "-m") {    
            p.max_m = checkInt(s, v);
        } else {
            cerr << s << " is not valid." << endl;
            exit(1);
        }
    }
    // check mandatory parameters
    if (p.input == "") {
        cerr << "-i is missing." << endl;
        exit(1);
    }
    if (p.output == "") {
        cerr << "-o is missing." << endl;
        exit(1);
    }
    // make sure output dir exits
    string cmd = "mkdir -p " + p.output;
    if(system(cmd.c_str())){}
}

int main(int argc, char** argv) {
    parameter p;
    checkParameter(argc, argv, p);
    omp_set_num_threads(p.thread);
    // run
    DetectIntersperse(p).run();
    return 0;
}
