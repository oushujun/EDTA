#ifndef _FASTAFILE_H_
#define _FASTAFILE_H_

#define CAPINC  64

typedef struct {
      int  id_len; /* allocated size of the sequence name string*/
      char    *id; /* id only, up to first space */
      int   d_len; /* allocated size of the comment string */
      char *descr; /* any comment on the defline, after the first space */
      int s_len; /* allocated length of the sequence string */
      int   len; /* the actual string length of seq */
      char* seq; /* the sequence itself */
} FastaSeq;

/* returns stdin if name is NULL and mode="r", 
   or stdout if name is NULL and mode is "w" */
FILE *myfopen(char* name, char* mode);

void *mymalloc(int size);
void *myrealloc(void* ptr, int size);
void check_eof(int c);
void bad_fastafmt(); /* invalid format encountered, exit */

/* allocates a FastaSeq structure on the heap and initializes it*/
FastaSeq* newFastaSeq(void);
void freeFastaSeq(FastaSeq* seq);

/* initializes a static FastaSeq structure */
void initFastaSeq(FastaSeq* seq);
/* clears (frees) internal strings of a FastaSeq structure */
void doneFastaSeq(FastaSeq* seq);


/* uses the given FastaSeq* seq storage structure; 
   can be called repeatedly for a multi-fasta file;
   returns NULL when no more records are found 
   othewise returns seq
  */
FastaSeq *getfasta(FILE* fp, FastaSeq* seq, int* is_last);

void putfasta(FastaSeq *fa, char* name);

#endif
