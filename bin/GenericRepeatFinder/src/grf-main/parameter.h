// Copyright (C) 2016 shij@miamioh.edu

#ifndef PARAMETER_H
#define PARAMETER_H

#include <string>
#include <climits>

using std::string;

struct parameter {
    int min_stem = 0;
    int min_tsd = 2;
    int max_tsd = 10;
    int match = 1;
    int mismatch = 1;
    int indel = 2;
    int percent = 10;
    int max_mismatch = INT_MAX;
    int max_indel = INT_MAX;
    int seed = 10;
    int thread = 1;
    int format = 0;
    int block = 100;
    int seed_mismatch = 1;
    int min_space = -1;
    int max_space = -1;
    int min_space_ir = 0;
    int min_space_mite = 30;
    int min_space_dr = 0;
    int max_space_ir = 980;
    int max_space_mite = 780;
    int max_space_dr = 980;
    int max_stem = INT_MAX;
    int min_spacer_len = 0;
    int max_spacer_len = INT_MAX;
    float r = -1;
    float br = 0.8;
    string input = "";
    string output = "";
    string choice = "";
    string program_path = "";
};

#endif /* PARAMETER_H */
