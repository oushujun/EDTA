/*
* stdaln.h -- standard alignment (local and banded global alignment)
*
* Copyright (c) 2003-2006, Heng Li <liheng@genomics.org.cn>
*                                  <lh3lh3@gmail.com>
*
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with this library; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*
*/

#ifndef LH3_STDALN_H_
#define LH3_STDALN_H_


#define STDALN_VERSION 0.9.4


#ifndef MYALLOC
# define MYALLOC malloc
#endif
#ifndef MYFREE
# define MYFREE free
#endif

#define FROM_M 0
#define FROM_I 1
#define FROM_D 2

/* This is the smallest integer. It might be CPU-dependent in very RARE cases. */
#define MINOR_INF -1073741823

typedef struct
{
    int gap_open;
    int gap_ext;
    int gap_end;

    int *matrix;
    int row;
    int band_width;
}
AlnParam;

typedef struct
{
    int i, j;
    unsigned char ctype;
}
path_t;

typedef struct
{
    path_t *path; /* for advanced users... :-) */
    int path_len; /* for advanced users... :-) */
    int start1, end1; /* start and end of the first sequence, coordinations are 1-based */
    int start2, end2; /* start and end of the second sequence, coordinations are 1-based */
    int score; /* score */

    char *out1, *out2; /* print them, and then you will know */
    char *outm;
}
AlnAln;

#ifdef __cplusplus
extern "C" {
#endif

    AlnAln *aln_stdaln_aux(const char *seq1, const char *seq2, const AlnParam *ap, int is_global, int len1, int len2);
    AlnAln *aln_stdaln(const char *seq1, const char *seq2, const AlnParam *ap, int is_global);
    void aln_free_AlnAln(AlnAln *aa);

#ifdef __cplusplus
}
#endif

/********************
 * global variables *
 ********************/

extern AlnParam aln_param_nt2nt; /* = { 10,  2,  2, aln_sm_nt, 16, 75 }; */
extern AlnParam aln_param_aa2aa; /* = { 20, 19, 19, aln_sm_read, 16, 75 }; */
extern AlnParam aln_param_rd2rd; /* = { 12,  2,  2, aln_sm_blosum62, 22, 50 }; */

/* common nucleotide score matrix for 16 bases */
extern int aln_sm_nt[];

/* BLOSUM62 and BLOSUM45 */
extern int aln_sm_blosum62[], aln_sm_blosum45[];

/* common read for 16 bases. note that read alignment is quite different from common nucleotide alignment */
extern int aln_sm_read[];

/* human-mouse score matrix for 4 bases */
extern int aln_sm_hs[];

#endif
