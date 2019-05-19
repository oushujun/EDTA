print ([[

File format for option '-offsetfile':

The file supplied to option '-offsetfile' defines a mapping table named
``offsets''. It maps the `sequence-region` entries given in the GFF3_file to
offsets.
It can be defined as follows:

    offsets = {
      chr1  = 1000,
      chr2  = 500
    }

When this example is used, all features with seqid ``chr1'' will be offset by
1000 and all features with seqid ``chr2'' by 500.

If '-offsetfile' is used, offsets for all `sequence-regions` contained in the
given GFF3 files must be defined.]])
