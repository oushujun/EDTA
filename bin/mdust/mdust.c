#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include "dust.h"
#include "fastafile.h"

/*
Based on the code of Tatusov & Lipman  (http://blast.wustl.edu/pub/dust/)
*/

#define usage "mdust [<fasta-file>] [-w <wsize>] [-v <cut-off>] [-m N|X|L] [-c] \n\
   if no <fasta-file> is given, a multi-fasta stream is expected at stdin\n\
   -v default <cut-off> value is 28 (lower values might mask more, \n\
      but possibly still useful sequence; > 64 will rarely mask poly-triplets)\n\
   -w set maximum word size to <wsize> (default 3)\n\
   -m if fasta output is not disabled by -c, set the masking letter type:\n\
        N ('N', default), X ('X'), L (make lowercase)\n\
   -c output masking coordinates only: \n\
        seq_name, seqlength, mask_start, mask_end  (tab delimited)\n\
"        

void err_param(char opt, char* value) {
 if (value!=NULL)
    fprintf(stderr, "%s\nInvalid option/value usage: '-%c %s'.\n", usage, opt, value);
   else
    fprintf(stderr, "%s\nInvalid option usage: '-%c'\n", usage, opt);
 exit(1);
}

int main(int argc, char* argv[]) {
      FastaSeq fa;
      DustRgn *rgn,*firstrgn;
      int wsize;
      FILE * inf;
      char* ffile=NULL;
      int level = 28; /* important: use 28 as a default value, to mask only long enough lc areas */
      int onlycoords=0;
      int i;
      char opt;
      char* value;
      int opst=1;
      maskFunc maskfn=&maskN;
      /*  parsing parameters: */
      if (argc>1) {
         if (argv[1][0]!='-'  ) {
            ffile=argv[1];
            opst++;
            }
          else if (argv[1][1]=='h' || argv[1][1]=='-' ||argv[1][1]=='\0') {
                 fprintf(stderr, "Usage: %s", usage);
                 return 1;
                 }
         }
       /* options only */  
       for (i=opst; i<argc;  i++) {
          if (argv[i][0]!='-') {
              fprintf(stderr, "Usage: %s", usage);
              exit(1);
              }
          opt=argv[i][1];    
          if (opt=='c') {
             onlycoords=1; 
             continue;
             }
          /* the rest of the options MUST have an associated value */
          if (strlen(argv[i])==2) {
              i++;
              if (i>=argc) err_param(opt, NULL);
              value=argv[i];
              }
             else {
              value=(argv[i]+2);
              }
          switch (opt) {
            case 'v': level=atoi(value);
                      if (level<=0) 
                         err_param(opt, value);
                      break;
            case 'w': wsize=atoi(value);
                      if (wsize<=0) 
                         err_param(opt, value);
                      set_dust_word(wsize);
                      break;
            case 'm': switch (tolower(*value)) {
                           case 'n':maskfn=&maskN;break;
                           case 'x':maskfn=&maskX;break;
                           case 'l':
                           case 'c':maskfn=&mask_lc;break;
                         default:
                          err_param(opt, value);
                         }
                      break;                      
            default: 
                err_param(opt,NULL);
                return 1;
             }
          }
          
      //set_dust_level(level);
      inf=myfopen(ffile,"r");
      initFastaSeq(&fa);
      
      if (onlycoords) {
           while (getfasta(inf,&fa, &i)!=NULL) {
               firstrgn = dust(&fa, level, NULL);
               for (rgn=firstrgn;rgn!=NULL; rgn=rgn->next) 
                  printf("%s\t%d\t%d\t%d\n", fa.id, fa.len, rgn->from+1, rgn->to+1);
               freeDustRgns(firstrgn);
               }
          }
       else 
         dust_file(inf, level, maskfn);
  doneFastaSeq(&fa);
  return 0;
}
