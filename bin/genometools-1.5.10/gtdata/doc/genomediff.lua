print ([[

The genomediff tool only accepts DNA input.

When used with sequence files or encseq, an enhanced suffix array will be
built in memory. The ESA will not be created completely, but construction will
use '-memlimit' as a threshold and build it partwise, calculating the Shu-length
for each part.

File format for option '-unitfile' (in Lua syntax):

  units = {
   genome1 = { "path/file1.fa", "file2.fa" },
   genome2 = { "file3.fa", "path/file4.fa" }
  }

Give the path to the files as they were given to the encseq tool!
You can use

  $ gt encseq info INDEXNAME

to get a list of files in an encoded sequence.

Comment lines in Lua start with '--' and will be ignored.

See `GTDIR/testdata/genomediff/unitfile1.lua` for an example.

Options '-pl', '-dc' and '-memlimit' are options to influence ESA construction.
]])
