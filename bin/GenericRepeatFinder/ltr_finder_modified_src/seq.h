/*
* nt_align, nucleotide-level alignment library
*
* Copyright (c) 2003-2004, Li Heng <liheng@genomics.org.cn>
*                                  <lihengsci@yahoo.com>
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

#ifndef SEQ_H_
#define SEQ_H_

#include <stdio.h>

#ifdef USE_KR_ALLOC
#include "../lib/mem.h"
#define MYALLOC kr_alloc
#define MYFREE kr_free
#define MYREALOC kr_realloc
#else
#include <malloc.h>
#define MYALLOC malloc
#define MYFREE free
#define MYREALLOC realloc
#endif /* USE_KR_ALLOC */


#define SEQ_BLOCK_SIZE 512
#define INIT_SEQ(seq) \
 (seq).s = 0; (seq).l = (seq).m = 0;
#define CHAR2QUAL(c) \
 ((isdigit(c))? ((c)-'0') : ((islower(c))? ((c)-'a'+10) : ((isupper(c))? ((c)-'A'+36) : 0)))
#define QUAL2CHAR(q) \
 (((q)<10)? ((q)+'0') : (((q)<36)? ((q)-10+'a') : (((q)<62)? ((q)-36+'A') : 'Z')))

typedef unsigned char uchar;

typedef struct
{
    int l, m; /* length and maximum buffer size */
    uchar *s; /* sequence */
}
seq_t;


int read_fasta(FILE*, seq_t*, char*, char*);
int read_qual(FILE*, seq_t*, char*, char*);

#endif /* SEQ_H_ */

