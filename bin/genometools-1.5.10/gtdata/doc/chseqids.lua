print ([[

File format for mapping_file:

The supplied mapping file defines a mapping table named ``chseqids''. It maps
the `sequence-region` entries given in the GFF3_file to other names.  It can be
defined as follows:

    chseqids = {
      chr1  = "seq1",
      chr2  = "seq2"
    }

When this example is used, all sequence ids ``chr1'' will be changed to ``seq1''
and all sequence ids ``chr2'' to ``seq2''.]])
