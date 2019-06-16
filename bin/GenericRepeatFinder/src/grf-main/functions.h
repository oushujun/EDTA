// Copyright (C) 2016 shij@miamioh.edu

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <string>
#include <vector>
#include <utility>

using std::string;
using std::vector;
using std::pair;

void fastaread(vector<string> & id, vector<string> & seq, const string & file);
double seqComplexity(const string & seq);
void mapVector(const string & seq, vector<pair<int, int> > & nseq);
void cumsum(vector<pair<int, int> > & v);
void readFile(vector<string> & v, const string & file);
int checkInt(const string & s, const string & v);
int checkInt2(const string & s, const string & v);
float checkFloat(const string & s, const string & v);
string compress(const string & s);
bool isPair(char a, char b);
string decompress(const string & s, char c);
void iniMatrix(vector<vector<pair<int, char > > > & matrix, int indel);
void iniMatrix2(vector<vector<pair<int, char > > > & matrix, int indel);
string reverseComplement(const string & s);
void fillMatrix(vector<vector<pair<int, char > > > & matrix,
        const string & a, const string & b, int match, int mismatch,
        int indel);
string reverse(const string & s);
int getStemLen(const string & s, char c);
string detectTSD(const string & seq, unsigned int start, unsigned int end, 
        int min, int max);
bool filterLowComplex(const string & seq, const string & tir);
float gcContent(const string & s);
void replace(string & s, char a);
void findError(vector<vector<pair<int, char > > > & matrix, int i, int j, 
        int & m2, int & i2, const string & a, const string & b);
bool match(char a, char b);

#endif /* FUNCTIONS_H */
