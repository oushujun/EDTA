// Copyright (C) 2016 shij@miamioh.edu

#ifndef PARAMETER_H
#define PARAMETER_H

using std::string;

struct parameter {
    unsigned thread = 1;
    unsigned format = 0;
    unsigned copy = 3;
    unsigned seed = 20;
    unsigned max_n = 1;
    unsigned max_m = 2;
    double min_con = 0.8;
    string input = "";
    string output = "";
};

#endif /* PARAMETER_H */
