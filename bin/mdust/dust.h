#ifndef _DUST_H_
#define _DUST_H_
#include "fastafile.h"

#define MAXREG    1001


typedef struct dustRgn{
     /* char* seqid; */
     int from;
     int to;
     struct dustRgn* next;
} DustRgn;


typedef char (*maskFunc)(char);

DustRgn* addDustRgn(DustRgn* last, int f, int t);
void freeDustRgns(DustRgn* first);

char maskN(char);
char maskX(char);
char mask_lc(char);

void set_dust_window(int value);
void set_dust_word(int value);

DustRgn* dust(FastaSeq* fa, int level, maskFunc maskfn);

void dust_file(FILE* f, int level, maskFunc maskfn);


#endif
