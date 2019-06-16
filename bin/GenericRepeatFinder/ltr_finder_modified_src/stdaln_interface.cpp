/*
* =====================================================================================
*
*       Filename:  stdaln_interface.cpp
*
*    Description:  interface to Stdaln lib
*
*        Version:  1.0
*        Created:  2006年11月21日 19时27分54秒 CST
*       Revision:  none
*       Compiler:  gcc
*
*         Author:   (), 
*        Company:  
*
* =====================================================================================
*/

#include "stdaln_interface.h"
#include "struct.h"
#include <iostream>

using namespace std;

int gap_open = 3;  //add by xz
int gap_ext = 1;
int gap_end = 1;
int score_match = 2;
int score_mismatch = -2;

int global_align(const char* str1, int begin1, int end1, const char* str2, int begin2, int end2, int* max_gap)
{
    AlnAln *aln_global;
    aln_global = aln_stdaln_aux(str1 + begin1, str2 + begin2, &aln_param_nt2nt, 1, end1 - begin1 + 1, end2 - begin2 + 1);
    int s = aln_global->score;
    int g_count = 0;
    int m_count = 0;
    *max_gap = 0;

    for (int i = 0;i < aln_global->path_len;++i)
    {
        if (aln_global->outm[i] == ' ') //miss match
        {
            g_count++;
        }
        else
        {
            if (g_count >= 2) //2 or more mismatch is gap, 1 is SNP//???here is unsure,may be we should check '-' in out1&2
            {

                if (CHECK_PAIRS)
                    cout << " gap between pairs:" << g_count << endl;

                *max_gap += g_count;
            }

            g_count = 0;
            ++m_count;
        }

    }
    aln_free_AlnAln(aln_global);
    //return s;
    return m_count;
}
/* initialize score matrix, only for nt2nt */
void set_score_matrix(int gap_open, int gap_ext, int gap_end, int match_score, int mismatch_score)
{
    int i, j, *p;
    int row = aln_param_nt2nt.row;

    for (i = 0; i < row; ++i)
    {
        for (j = 0; j < row; ++j)
        {
            p = aln_param_nt2nt.matrix + i * row + j;

            if (*p == 2)
            { /* match */
                *p = score_match; //by xz
            }
            else if (*p == -2)
            {
                *p = score_mismatch; //by xz, else *p = -2;
            }
            //printf("%4d",*p);
        }
        //printf("\n");
    }

    aln_param_nt2nt.gap_open = gap_open;
    aln_param_nt2nt.gap_ext = gap_ext;
    aln_param_nt2nt.gap_end = gap_end;
}

