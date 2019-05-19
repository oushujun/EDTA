print ([[

If neither option '-check' nor option '-duplicates' is used, the fingerprints
for all sequences are shown on stdout.

Fingerprint of a sequence is case insensitive. Thus MD5 fingerprint of two
identical sequences will be the same even if one is soft-masked.

Examples
--------

Compute (unified) list of fingerprints:

    $ gt fingerprint U89959_ests.fas | sort | uniq > U89959_ests.checklist_uniq

Compare fingerprints:

    $ gt fingerprint -check U89959_ests.checklist_uniq U89959_ests.fas
    950b7715ab6cc030a8c810a0dba2dd33 only in sequence_file(s)

Make sure a sequence file contains no duplicates (not the case here):

    $ gt fingerprint -duplicates U89959_ests.fas
    950b7715ab6cc030a8c810a0dba2dd33        2
    gt fingerprint: error: duplicates found: 1 out of 200 (0.500%)

Extract sequence with given fingerprint:

    $ gt fingerprint -extract 6d3b4b9db4531cda588528f2c69c0a57 U89959_ests.fas
    >SQ;8720010
    TTTTTTTTTTTTTTTTTCCTGACAAAACCCCAAGACTCAATTTAATCAATCCTCAAATTTACATGATAC
    CAACGTAATGGGAGCTTAAAAATA

Return values
-------------

- 0  everything went fine ('-check': the comparison was successful;
                           '-duplicates': no duplicates found)
- 1  an error occurred     ('-check': the comparison was not successful;
                           '-duplicates': duplicates found)]])
