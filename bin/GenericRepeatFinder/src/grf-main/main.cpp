// Copyright (C) 2016 shij@miamioh.edu

#include <cstdlib>
#include <iostream>
#include <string>
#include <exception>
#include <ctime>
#include <omp.h>
#include "parameter.h"
#include "functions.h"
#include "DetectIR.h"
#include "DetectDR.h"
#include "DetectMITE.h"

using namespace std;

void checkParameter(int argc, char** argv, parameter & p) {
    string name = "grf-main";
    string help =
            "This program performs genome-wide terminal inverted repeat (TIR), terminal direct repeat (TDR), "
            "and MITE candidate detection for input genome sequences.\n"
            "Usage: " + name + " [options]\n"
            "[mandatory parameters]\n"
            "-i <string>  Input genome sequence file in FASTA format.\n"
            "-c <int>  Choice of analysis: 0: TIR detection; 1: MITE candidate detection; 2: TDR detection.\n"
            "-o <string>  Output directory.\n"
            "--min_tr <int>  Minimum length of the terminal repeats; must >= 5.\n"
            "[optional parameters]\n"
            "-t <int>  Number of threads used in this program; default = 1.\n"
            "-f <int> Format for outputs; 0: FASTA format; 1: only IDs; default = 0.\n"
            "--max_mismatch <int>  Maximum number of mismatches allowed in the terminal repeats; "
            "set -1 to be unlimited; default = -1.\n"
            "--max_indel <int>  Maximum number of indels allowed in the terminal repeats; "
            "set -1 to be unlimited; default = -1.\n"
            "-p <int>  Maximum percentage of unpaired nucleotides in the terminal repeats; "
            "set -1 to be unlimited; default = 10.\n"
			"-r <float>  Maximum length ratio of spacer/total sequence; set -1 to be unlimited; "
			"default = -1.\n"
            "-s <int>  Length of the seed region; default = 10; must >= 5 and <= '--min_tr'.\n"
            "--seed_mismatch <int>  Maximum mismatch number in the seed region; default = 1.\n"
            "--min_space <int>  Minimum distance between two seed regions; "
            "for TIRs/TDRs, default = 0; for MITEs, default = 30.\n"
            "--max_space <int>  Maximum distance between two seed regions; "
            "for TIRs/TDRs, default = 980; for MITEs, default = 780.\n"
            "\n"
            "If indel is enabled, in the extension of seed regions, the following scoring matrix for alignment is used.\n"
            "--match <int>  Award score (positive number) for 1 match; default = 1.\n"
            "--mismatch <int>  Penalty score (positive number) for 1 mismatch; default = 1.\n"
            "--indel <int>  Penalty score (positive number) for 1 indel; default = 2.\n"
            "--block <int>  Block size during alignment; default = 100.\n"
            "--block_ratio <float>  For the best alignment in the current block, "
            "if the length of aligned sequences <= block_ratio * block_size, the alignment procedure will stop "
            "and the end position of the best alignment will be returned. "
            "Otherwise, a new block will be created and the alignment will continue from the current end position; "
            "default = 0.8.\n"
            "\n"
            "For MITE detection,\n"
            "--min_tsd <int>  Minimum length of TSDs; default = 2.\n"
            "--max_tsd <int>  Maximum length of TSDs; default = 10.\n"
            "\n"
            "To restrict terminal repeat (TR) and spacer length of detected TIRs, TDRs, and MITE candidates, set the following options:\n"
            "Note: minimum TR length has been set in the mandatory parameter '--min_tr'.\n"
            "--max_tr <int>  Maximum TR length.\n"
            "--min_spacer_len <int>  Minimum spacer length.\n"
            "--max_spacer_len <int>  Maximum spacer length.\n"
            ;
    if (argc == 1 || string(argv[1]) == "-h") {
        cout << help << endl;
        exit(0);
    }
    
    // get program path
    {
        string s = argv[0];
        size_t pos = s.find_last_of("/");
        if (pos == string::npos) {
            p.program_path = ".";
        } else {
            p.program_path = s.substr(0, pos);
        }
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
        } else if (s == "-c") {
            p.choice = v;
        } else if (s == "-o") {
            p.output = v;
        } else if (s == "-r") {
            p.r = checkFloat(s, v);
            if (p.r < 0 && p.r != -1) {
                cerr << s << " " << v << " is invalid." << endl;
                exit(1);
            }
        } else if (s == "--min_tr") {
            p.min_stem = checkInt(s, v);
            if (p.min_stem < 5) {
                cerr << s << " " << v << " is out of range." << endl;
                exit(1);
            }
        } else if (s == "--max_mismatch") {
            if (v != "-1") {
                p.max_mismatch = checkInt2(s, v);
            }
        } else if (s == "--max_indel") {
            if (v != "-1") {
                p.max_indel = checkInt2(s, v);
            }
        } else if (s == "-p") {
            if (v == "-1") {
                p.percent = 100;
            } else {
                p.percent = checkInt2(s, v);
                if (p.percent >= 40) {
                    cerr << s << " " << v << " is too high." << endl;
                    exit(1);
                }
            }
        } else if (s == "--match") {
            p.match = checkInt2(s, v);
        } else if (s == "--mismatch") {
            p.mismatch = checkInt2(s, v);
        } else if (s == "--indel") {            
            p.indel = checkInt2(s, v);
        } else if (s == "--block") {
            p.block = checkInt(s, v);
        } else if (s == "--block_ratio") {
            p.br = checkFloat(s, v);
            if (p.br <= 0 || p.br > 1) {
                cerr << s << " " << v << " is invalid." << endl;
                exit(1);
            }
        } else if (s == "-s") {
            p.seed = checkInt(s, v);
            if (p.seed < 5) {
                cerr << s << " " << v << " is out of range." << endl;
                exit(1);
            }
        } else if (s == "--seed_mismatch") {
            p.seed_mismatch = checkInt2(s, v);
        } else if (s == "--min_space") {
            p.min_space = checkInt2(s, v);
        } else if (s == "--max_space") {
            p.max_space = checkInt2(s, v); 
        } else if (s == "--min_tsd") {
            p.min_tsd = checkInt(s, v);
        } else if (s == "--max_tsd") {
            p.max_tsd = checkInt(s, v);
        } else if (s == "-t") {
            p.thread = checkInt(s, v);
        } else if (s == "-f") {
            p.format = checkInt2(s, v);
            if (p.format != 0 && p.format != 1) {
                cerr << s << " " << v << "is not valid." << endl;
                exit(1);
            }        
        } else if (s == "--max_tr") {
            p.max_stem = checkInt2(s, v);
        } else if (s == "--min_spacer_len") {
            p.min_spacer_len = checkInt2(s, v);
        } else if (s == "--max_spacer_len") {
            p.max_spacer_len = checkInt2(s, v);
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
    if (p.choice == "") {
        cerr << "-c is missing." << endl;
        exit(1);
    }
    if (p.output == "") {
        cerr << "-o is missing." << endl;
        exit(1);
    }
    if (p.min_stem == 0) {
        cerr << "--min_tr is missing." << endl;
        exit(1);
    }
    // check seed len and mismatch
    if (p.seed > p.min_stem) {
        cerr << "'-s' should be <= '--min_tr'." << endl;
        exit(1);
    }
    if (p.seed_mismatch > p.max_mismatch) {
        cerr << "'--seed_mismatch' should be <= '--max_mismatch'." << endl;
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
    if (p.choice == "0") {
        // inverted repeat
        if (p.min_space == -1) {
            p.min_space = p.min_space_ir;
        }
        if (p.max_space == -1) {
            p.max_space = p.max_space_ir;
        }
        DetectIR(p).run();
    } else if (p.choice == "1") {
        // mite
        if (p.min_space == -1) {
            p.min_space = p.min_space_mite;
        }
        if (p.max_space == -1) {
            p.max_space = p.max_space_mite;
        }
        DetectMITE(p).run();
    } else if (p.choice == "2") {
        // direct repeat
        if (p.min_space == -1) {
            p.min_space = p.min_space_dr;
        }
        if (p.max_space == -1) {
            p.max_space = p.max_space_dr;
        }
        DetectDR(p).run();
    } else {
        cerr << "-c is not valid." << endl;
        exit(1);
    }
    // end
    time(&rawtime);
    cout << "end: " << ctime(&rawtime) << endl;
    return 0;
}
