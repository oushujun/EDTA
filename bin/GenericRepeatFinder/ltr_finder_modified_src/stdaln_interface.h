/*
* =====================================================================================
* 
*       Filename:  stdaln_interface.h
* 
*    Description:  interface to Stdaln lib
* 
*        Version:  1.0
*        Created:  2006年11月21日 18时50分22秒 CST
*       Revision:  none
*       Compiler:  gcc
* 
*         Author:  XuZhao (), 
*        Company:  
* 
* =====================================================================================
*/
#ifndef STDALN_INTERFACE_H
#define STDALN_INTERFACE_H

#include "stdaln.h"
extern int gap_open ; // 6;
extern int gap_ext ; // 1;
extern int gap_end ; // 2;
extern int score_match ; // 2;
extern int score_mismatch ; // -5;

int global_align(const char* str1, int begin1, int end1,
                 const char* str2, int begin2, int end2, int* max_gap);
void set_score_matrix(int gap_open, int gap_ext, int gap_end,
                      int match_score, int mismatch_score);

#endif
