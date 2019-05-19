print ([[

Example:
--------

Let's assume we have a GFF3 file 'csa_example_spliced_alignments.gff3'
containing the following four overlapping spliced alignments (represented as
genes with exons as children):
]])
print(io.open(gtdata_doc_dir.."csa_example_spliced_alignments.gff3"):read("*a"))
print([[
To compute the consensus spliced alignments we call:

    $ gt csa csa_example_spliced_alignments.gff3

Which returns:
]])
print(io.open(gtdata_doc_dir.."csa_example_consensus_spliced_alignments.gff3"):read("*a"))
print([[
As one can see, they have been combined into a consensus spliced alignment
(represented as genes with mRNAs as children which in turn have exons as
children) with two alternative splice forms. The first and the third spliced
alignment have been combined into the first alternative splice form (mRNA1) and
the the second and the fourth spliced alignment into the second alternative
splice form (mRNA2).

As one can see, the second exon from the first alternative splice form is
shorter than the corresponding exon from the second alternative splice form.]])
