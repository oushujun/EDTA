#include <stdio.h>
#include <ctype.h>
#include <malloc.h>

#include "fastafile.h"


#define TMPNAME  "tmp.seg"

void myabort(char* mess) {
      fprintf(stderr, "Error: %s.\n", mess);
      exit(1);
}

void *mymalloc(int size) {
      void *buf;
      if ((buf = malloc(size)) == NULL) {
            myabort("Not enough memory");
      }
      return buf;
}

void *myrealloc(void* buf, int size) {
      void *result;
      if ((result = realloc(buf, size)) == NULL) {
            myabort("Not enough memory");
      }
      return result;
}

/* allocates a FastaSeq structure on the heap and initializes it*/
FastaSeq* newFastaSeq(void) {
 FastaSeq* r=(FastaSeq*)mymalloc(sizeof(FastaSeq));
 initFastaSeq(r);
 return r;
 }

void freeFastaSeq(FastaSeq* seq) {
 doneFastaSeq(seq);
 free(seq);
 }

/* initializes a static FastaSeq structure */
void initFastaSeq(FastaSeq* seq) {
 seq->id=(char*)mymalloc(CAPINC);
 seq->id_len=CAPINC;
 seq->id[0]='\0';
 seq->descr=(char*)mymalloc(CAPINC);
 seq->descr[0]='\0';
 seq->d_len=CAPINC;
 seq->seq=(char*)mymalloc(CAPINC<<1);
 seq->seq[0]='\0';
 seq->len=0;
 seq->s_len=CAPINC<<1; 
 }
 
/* clears (frees) internal strings of a FastaSeq structure */
void doneFastaSeq(FastaSeq* seq) {
 free(seq->id); seq->id=NULL; seq->id_len=0;
 free(seq->descr); seq->descr=NULL; seq->d_len=0;
 free(seq->seq); seq->seq=NULL; seq->s_len=0;seq->len=0;
 }

 

/* initializes a static FastaSeq structure */
void initFastaSeq(FastaSeq* seq);
/* clears (frees) internal strings of a FastaSeq structure */
void doneFastaSeq(FastaSeq* seq);


FILE *myfopen(char* name, char* mode) {
      FILE *fp;

      if (name == NULL) {
            return (mode != NULL && *mode == 'r') ? stdin : stdout;
      }
      if ((fp = fopen(name, mode)) == NULL) {
            if (mode != NULL && *mode == 'r') {
                  fprintf(stderr, "No such file: \"%s\"\n", name);
            } else {
                  fprintf(stderr, "Failed to open file: \"%s\"\n", name);
            }
            exit(1);
      }
      return fp;
}

void bad_fastafmt() {
      myabort("Not a FastaSeq file");
}

void check_eof(int c) {
      if (c == EOF) 
            bad_fastafmt();
}


/* the first character must be '>' for each call
   seq must be an initialized FastaSeq structure
   */
  
FastaSeq *getfasta(FILE* fp, FastaSeq* seq, int* is_last) {
      int c, len; 
      int* buflen; 
      char** buf;
      int before;      
      /* fp = myfopen(name, "r"); */
      c = getc(fp);
      if (c==EOF) return NULL;
      if (c != '>')
            bad_fastafmt();
      len = 0; //chars accumulated so far
      /* -------- read the defline first */
      buflen=&seq->id_len;
      buf=&seq->id;
      before=1;
      while ((c = getc(fp)) != EOF && c != '\n') {
          if (len >= *buflen-1) {
                  *buf=(char*)myrealloc(*buf, *buflen + CAPINC);
                  *buflen+=CAPINC;
                  }
          if (before && isspace(c)) {
             /* space encountered => seq_name finished */
             before=0;
             (*buf)[len]='\0';
             buf=&seq->descr;
             buflen=&seq->d_len;             
             len=0;
             continue; /* skip this space */
             }
          (*buf)[len]=c;
          len++;
          }
      (*buf)[len]='\0'; /* terminate the comment string */    
      check_eof(c); /* it's wrong to have eof here */
      /*----- read the sequence now: */
      len=0;
      before=1; //newline before indicator
      while ((c = getc(fp)) != EOF && c != '>') {
            if (isspace(c)) {
                   before = (c=='\n')?1:0;
                   continue; /* skip spaces */
                   }
            if (len >= seq->s_len-1) {
                  seq->seq = (char*)myrealloc(seq->seq, seq->s_len + CAPINC);
                  seq->s_len+=CAPINC;
                  }
            seq->seq[len] = c;
            before=0;
            len++;
            }
      seq->seq[len] = '\0';
      seq->len=len;      
      if (c=='>') {
         if (!before) bad_fastafmt(); /* '>' must only be at start of line, 
                                       never within the sequence ! */
         *is_last=0; /* FALSE - not the last one */
         ungetc(c, fp);
         }
        else  *is_last=1; /* TRUE - eof() here */
      return seq;
}

void putfasta(FastaSeq *fa, char* name) {
      FILE *fp;
      char *s;
      int i, i60;

      fp = myfopen(name, "w");
      s = (fa->id == NULL) ? (char*)"ANONYMOUS" : fa->id;
      if (*s != '>') 
            putc('>', fp);
      fwrite(s, 1, strlen(s), fp);
      i=strlen(fa->descr);
      if (i>0) {
        putc(' ',fp);
        fwrite(fa->descr, 1, i, fp);
        }
      i60 = 60;      
      for (i=0, s = fa->seq; i < fa->len; i++, s++, i60++) {
            if (i60 == 60) {
                  putc('\n', fp);
                  i60 = 0;
                  }
            putc(*s, fp);
      }
      putc('\n', fp);
      if (fp != stdout) {
            fclose(fp);
      }
}


