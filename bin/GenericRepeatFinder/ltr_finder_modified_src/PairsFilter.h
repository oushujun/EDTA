/*
* =====================================================================================
* 
*       Filename:  PairsFilter.h
* 
*    Description:  
* 
*        Version:  1.0
*        Created:  2006年10月20日 09时52分36秒 CST
*       Revision:  none
*       Compiler:  gcc
* 
*         Author:   (), 
*        Company:  
* 
* =====================================================================================
*/
#ifndef PAIRFILTER_H_LTR_FINDER
#define  PAIRFILTER_H_LTR_FINDER

#include "struct.h"
#include "seq.h"
#include "PBS.h"
#include <set>
#include <vector>

using namespace std;
int GetPairs (char *myString, int stringLength, int Lex, int Lmax, int Dmin, int Dmax,
              vector < candidate > &Pair, bool LeftMaxAlign = true);
void JoinPairs(char* myString, int stringLength, vector<candidate>& Pair, vector<stick>& sticks);
void EnumPossibleLTR (char* myString, int stringLength, PBS& pbs, PPT& ppt,
                      stick & st, vector < stick > &pLTRs);
bool ExtendPairs(char* myString, int stringLength, stick& st);
void FindSignal(char* myString, int maxLen, PBS& pbs,
                PPT& ppt, stick& st);
int FindEdge(char* str, int len, int win_size, int direct);
int FindTSR(char* myString, int pos1, int pos2);
void FindTwoChar(AlnAln *aln, vector<int>& pos1,
                 vector<int>& pos2, char* sub);
void ScoreCompensate(char* myString, TG_CA_TSR& tct);
void outstr(char *str, int pos, int len);
void EraseOverlap(vector<stick>& sticks);
void ShorterAlignStr(string& str);
void RefineSimilarity(char* str, char* aux, int len, int pos1, int pos2, stick& st);
#endif
