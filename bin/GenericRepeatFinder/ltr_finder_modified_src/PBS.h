/*
* =====================================================================================
* 
*       Filename:  PBS.h
* 
*    Description:  class for finding PBS
* 
*        Version:  1.0
*        Created:  2006年10月22日 22时58分23秒 CST
*       Revision:  none
*       Compiler:  gcc
* 
*         Author:   (), 
*        Company:  
* 
* =====================================================================================
*/

#ifndef PBS_H_LTR_FINDER
#define PBS_H_LTR_FINDER

#include <string>
#include <vector>

using namespace std;

class PBS
{
    vector< string > tRNA;
    vector< string > _tRNA;
    vector< string > tRNA_name;
    vector< string > _tRNA_name;

    int PBS_len;
    int minLen;

    void init();

public:
    PBS();
    PBS(char* filename);
    ~PBS(){
    };
    string LoadSeq(char* filename);
    int Search(const char *str, int len, string& name, int& begin, int& end, int& match_len, string& outstr);
    int rSearch(const char *str, int len, string& name, int& begin, int& end, int& match_len, string& outstr);
};

class PPT
{
    int window_size;
    int delta;
    void init();


public:
    PPT();
    PPT(int w, int d);
    ~PPT(){
    };
    int Search(char *str, int len, int& begin, int& match_count);
    int rSearch(char *str, int len, int& begin, int& match_count);
};
#endif
