print ([[

Correction mode:

One of the options '-aggressive', '-moderate', '-conservative' or '-expert'
must be selected.

The '-aggressive', '-moderate' and '-conservative' modes are presets of
the criteria by which it is decided if an observed discrepancy in
homopolymer length between cognate sequence and a read shall be corrected
or not. A description of the single criteria is provided by using
the '-help+'' option. The presets are equivalent to the following settings:

                        -aggressive    -moderate      -conservative
    -hmin               3              3              3
    -read-hmin          1              1              2
    -altmax             1.00           0.99           0.80
    -refmin             0.00           0.00           0.10
    -mapqmin            0              10             21
    -covmin             1              1              1
    -clenmax            unlimited      unlimited      unlimited
    -allow-multiple     yes            yes            no

The aggressive mode tries to maximize the sensitivity, the conservative
mode to minimize the false positives. An even more conservative set
of corrections can be achieved using the '-ann' option (see '-help+').

The '-expert' mode allows one to manually set each parameter; the default
values are the same as in the '-conservative' mode.

(Finally, for evaluation purposes only, the '-state-of-truth' mode can be used:
this mode assumes that the sequenced genome has been specified
as cognate sequence and outputs an ideal list of corrections.)]])
