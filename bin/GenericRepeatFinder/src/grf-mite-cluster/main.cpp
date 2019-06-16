// Copyright (C) 2016 shij@miamioh.edu

#include <cstdlib>
#include <iostream>
#include <string>
#include <omp.h>
#include <ctime>
#include "FamilyCluster.h"

using std::cout;
using std::cerr;
using std::endl;

int checkInt(const string & s, const string & v) {
    int i;
    try {
        i = stoi(v);
    } catch (std::exception & e) {
        cerr << s << " is invalid." << endl;
        exit(1);
    }
    if (i <= 0) {
        cerr << s << " is invalid." << endl;
        exit(1);
    }
    return i;
}

int checkInt2(const string & s, const string & v) {
    int i;
    try {
        i = stoi(v);
    } catch (std::exception & e) {
        cerr << s << " is invalid." << endl;
        exit(1);
    }
    if (i < 0) {
        cerr << s << " is invalid." << endl;
        exit(1);
    }
    return i;
}

void checkParameter(int argc, char** argv, parameter & p) {
    string name = "grf-mite-cluster";
    string help =
            "This program performs family clustering and filtration for MITE candidates.\n"
            "Usage: " + name + " [options]\n"
            "[mandatory parameters]\n"
            "-i <string>  'cd-hit-est' output file: \"*.clstr\".\n"
            "-g <string>  Input genome sequence file in FASTA format.\n"
            "-o <string>  Output directory.\n"
            "\n"
            "[optional parameters]\n"
            "-t <int>  Number of threads used in this program; default = 1.\n"
            "-f <int>  Length of flanking sequences of MTIEs; default = 50.\n"
            "-c <int>  Minimum copy number of MITEs in the genome; default = 3.\n"
            "In comparing flanking sequences of MITE candidates, the scoring matrix for alignment is as below:\n"
            "--match <int>  Award score (positive number) for 1 match; default = 1.\n"
            "--mismatch <int>  Penalty score (positive number) for 1 mismatch; default = 1.\n"
            "--indel <int>  Penalty score (positive number) for 1 indel; default  = 2.\n"
            ;
    if (argc == 1 || string(argv[1]) == "-h") {
        cout << help << endl;
        exit(0);
    }
    for (int i = 1; i < argc; i += 2) {
        string s(argv[i]);
        if (i + 1 >= argc) {
            cerr << "total parameter number not correct" << endl;
            exit(1);
        }
        string v(argv[i + 1]);
        if (s == "-i") {
            p.input = v;
        } else if (s == "-g") {
            p.genome = v;
        } else if (s == "-o") {
            p.output = v;
        } else if (s == "--match") {
            p.match = checkInt2(s, v);
        } else if (s == "--mismatch") {
            p.mismatch = checkInt2(s, v);
        } else if (s == "--indel") {
            p.indel = checkInt2(s, v);
        } else if (s == "-f") {
            p.flank = checkInt(s, v);
        } else if (s == "-t") {
            p.thread = checkInt(s, v);
        } else if (s == "-c") {
            p.copy = checkInt(s, v);    
        } else {
            cerr << s << " is not valid." << endl;
            exit(1);
        }
    }
    if (p.input == "") {
        cerr << "-i <string> is missing." << endl;
        exit(1);
    }
    if (p.genome == "") {
        cerr << "-g <string> is missing." << endl;
        exit(1);
    }
    if (p.output == "") {
        cerr << "-o <string> is missing." << endl;
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
    // start
    time_t rawtime;
    time(&rawtime);
    cout << "start: " << ctime(&rawtime) << endl;
    FamilyCluster(p).run();
    // end
    time(&rawtime);
    cout << "end: " << ctime(&rawtime) << endl;
    return 0;
}
