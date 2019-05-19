print ([[

The option -keys allows one to extract substrings or sequences from the given
sequence file or from a fasta index.
The substrings to be extracted are specified in a key file given
as argument to this option. The key file must contain lines of the form

  k

or

  k i j

where k is a string (the key) and the optional i and j are positive integers
such that i<=j. k is the key and the optional numbers i and j specify the
first position of the substring and the last position of the substring to be
extracted. The positions are counted from 1. If k is identical to the string
between the first first and second occurrence of the symbol | in a fasta
header, then the fasta header and the corresponding sequence is output.
For example in the fasta header

  >tr|A0AQI4|A0AQI4_9ARCH Putative ammonia monooxygenase (Fragment)

the fasta key is A0AQI4. If i and j are both specified, then the corresponding
substring is shown in fasta format. In the latter case the header of the
fasta formatted sequence in the output begins with

  >k i j

followed by the original original fasta header.

If the sequence input are fasta files, then the following holds:

  - duplicated lines in the input file lead to only one sequence in the output
  - the sequences are output according to the order in the original sequence
    files
  - the formatting of the output can be controlled by the options '-width',
    '-o', '-gzip', and '-bzip2'

If the sequence input comes from a fasta index (see below), the following holds:

  - option '-width' is required
  - option '-o', '-gzip' and '-bzip2' do not work
  - the sequences are output in the order the corresponding keys appear in
    the key file

If the end of the argument list only contains one filename, say fastaindex, then
it is checked if there is a file `fastaindex.kys`. This makes up part of the
fasta index, which is constructed by calling the suffixerator tool as follows:

  gt suffixerator -protein -ssp -tis -des -sds -kys -indexname fastaindex \
    -db inputfile1 [inputfile2 ..]

This reads the protein sequence files given to the option '-db' and creates
several files:

 - a file `fastaindex.esq` representing the sequence.
 - a file `fastaindex.ssp` specifying the sequence separator positions.
 - a file `fastaindex.des` showing the fasta headers line by line.
 - a file `fastaindex.sds` giving the sequence header delimiter positions.
 - a file `fastaindex.kys` containing the keys in the fasta files.

For the suffixerator command to work, the keys of the form |key| in the fasta
header must satisfy the following constraints:

  - they all have to be of the same length, not longer than 128, and not shorter
    than 1
  - they have to appear in lexicographic order]])
