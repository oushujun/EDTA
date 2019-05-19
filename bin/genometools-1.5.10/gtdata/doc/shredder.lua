print ([[

Each sequence given in 'sequence_file' is shreddered into consecutive pieces of
random length (between '-minlength' and '-maxlength') until it is consumed.
By this means the last shreddered fragment of a given sequence can be shorter
than the argument to option '-minlength'.
To get rid of such fragments use `gt seqfilter` (see example below).

Examples:
---------

Shredder a given BAC:

    $ gt shredder U89959_genomic.fas > fragments.fas

Shredder an EST collection into pieces between 50 and 100 bp and get rid of all
(terminal) fragments shorter than 50 bp:

    $ gt shredder -minlength 50 -maxlength 100 U89959_ests.fas \
      | gt seqfilter -minlength 50 - > fragments.fas
    # 130 out of 1260 sequences have been removed (10.317%)

Shredder an EST collection and show only random 10% of the resulting fragments:

    $ gt shredder -sample 0.1 U89959_ests.fas]])
