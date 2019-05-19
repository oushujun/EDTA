#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include "dust.h"
#include "fastafile.h"

#define TRUE   1
#define FALSE  0


static int word = 3; 
static int window = 64; 
static int window2 = 32; 

static int mv, iv, jv;

void set_dust_window(int value) {
      window = value;
      window2 = window / 2;
}

void set_dust_word(int value) {
      word = value;
}


DustRgn* addDustRgn(DustRgn* prev, int f, int t) {  
  DustRgn* r;
  int nomerge;  
  nomerge=1;
  /* fprintf(stderr, "addDustRgn(, %d, %d)\n", f, t); */
  if (t<0 || f>=t) return prev; /* not a maskable range */
  if (prev!=NULL) {
    if (f>prev->from && f<=prev->to) { //extend f
          if (prev->to<t) prev->to=t;
          nomerge=0;
          }
    if (prev->from>=f && prev->from<=t) {
          if (f<prev->from) prev->from=f;
          nomerge=0;
          }
    } /* prev!=NULL */    
  if (nomerge) { /* no merge, new range */
      r=(DustRgn*)malloc(sizeof(DustRgn));
      if (prev!=NULL) prev->next=r;      
      r->from=f;r->to=t;
      r->next=NULL;
      }
     else r=prev;    
  return r;
  }

void freeDustRgns(DustRgn* first) {
  DustRgn* p=first;
  DustRgn* next;
  while (p!=NULL) {
    next=p->next;
    free(p);
    p=next;
    }
  }


static void wo1(int len, char* s, int ivv) {
      int i, ii, j, v, t, n, n1, sum;
      static int counts[32*32*32];
      static int iis[32*32*32];
      int js, nis;

      n = 32 * 32 * 32;
      n1 = n - 1;
      nis = 0;
      i = 0;
      ii = 0;
      sum = 0;
      v = 0;
      for (j=0; j < len; j++, s++) {
            ii <<= 5;
            if (isalpha(*s)) {
                  if (islower(*s)) {
                        ii |= *s - 'a';
                  } else {
                        ii |= *s - 'A';
                  }
            } else {
                  i = 0;
                  continue;
            }
            ii &= n1;
            i++;
            if (i >= word) {
                  for (js=0; js < nis && iis[js] != ii; js++) ;
                  if (js == nis) {
                        iis[nis] = ii;
                        counts[ii] = 0;
                        nis++;
                  }
                  if ((t = counts[ii]) > 0) {
                        sum += t;
                        v = 10 * sum / j;
                        if (mv < v) {
                              mv = v;
                              iv = ivv;
                              jv = j;
                        }
                  }
                  counts[ii]++;
            }
      }
}

static int wo(int len, char* s, int* beg, int* end) {
      int i, l1;

      l1 = len - word + 1;
      if (l1 < 0) {
            *beg = 0;
            *end = len - 1;
            return 0;
      }
      mv = 0;
      iv = 0;
      jv = 0;
      for (i=0; i < l1; i++) {
            wo1(len-i, s+i, i);
      }
      *beg = iv;
      *end = iv + jv;
      return mv;
}

char maskN(char c) {
 return 'N';
}

char maskX(char c) {
 return 'X';
}

char mask_lc(char c) {
 return (char)tolower(c);
}


DustRgn* dust(FastaSeq* fa, int level, maskFunc maskfn) {
      int i, j, l, from, to, a, b, v;
      DustRgn* first=NULL;
      DustRgn* cur=NULL;
      DustRgn* prev=NULL;
      from = 0;
      to = -1;
      a=0;b=0;
      for (i=0; i < fa->len; i += window2) {
            from -= window2;
            to -= window2;
            l = (fa->len > i+window) ? window : fa->len-i;
            v = wo(l, fa->seq+i, &a, &b);
            if (maskfn==NULL) { /* coordinates only */
               /* return coordinates from-to*/
               /*fprintf(stderr, "%s : %d - %d \n", fa->id, i+from, i+to); */
               prev=cur;
               cur=addDustRgn(cur, i+from, i+to);
               if (first==NULL && cur!=prev) first=cur;
               }
             else {
               for (j = from; j <= to; j++) {
                     if (isalpha(fa->seq[i+j]))
                         fa->seq[i+j] = (*maskfn)(fa->seq[i+j]);
                     }
               }
            if (v > level) {
               if (maskfn==NULL) { 
                   /* return coordinates from a to min(b,window2)*/
                  /* fprintf(stderr, "%s : %d - %d \n", fa->id, 
                        i+a, i+((window2>b)? b : window2)); */
                  prev=cur;      
                  cur=addDustRgn(cur, i+a, i + ((window2>b)? b : window2));
                  if (first==NULL && cur!=prev) first=cur;
                  }
                else 
                  for (j = a; j <= b && j < window2; j++) {
                        if (isalpha(fa->seq[i+j]))
                              fa->seq[i+j] = (*maskfn)(fa->seq[i+j]);
                        }
                from = (window2>b)? b : window2;
                to = b;
                }
             else {
                from = 0;
                to = -1;
                }
      }
 return first;     
}

void dust_file(FILE *fin, int level, maskFunc maskfn) {
  int is_last;
  FastaSeq seq;
  initFastaSeq(&seq);
  while(getfasta(fin,&seq, &is_last)!=NULL) {
    dust(&seq, level, maskfn);
    putfasta(&seq, NULL);
    }
  doneFastaSeq(&seq);
}

