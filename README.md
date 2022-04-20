Fqzcomp5
========

This is a reimplementation of the old fqzcomp-4.x tool using the new
htscodecs library (with a few amendments, to be folded back upstream).

Given the parameterisation of the fqzcomp qual model, this is
currently a little slower than the v4.6 code, but this may potentially
be fixable by compiling custom implementations for regularly used
configurations.

The file format is also now block based, typically between 100MB and
1GB in size.  This permits it to be more aggressively multi-threaded,
but it harms the maximum compression ratio slightly.  In principle
this also permits some level of random access, which for fastq
practically means sub-dividing the data into blocks, e.g. for
parallel dispatching to an aligner.  As yet the index has not been
written to permit this.

Compared to fqzcomp-4:

- New fast modes (-1 and -3) using basic bit-packing and rANS.  This
  is designed for rapid fastq compression, and is also appropriate for
  smaller block sizes.

- More flexible fqzcomp_qual modelling, which on some data sets can
  give substantially smaller quality compression depending on the
  nature of the data.

- The optional addition of sequence based in the quality compression
  context.  This improves compression ratios for ONT and PacBio data
  sets.

- A learning algorithm (similar to our CRAM implementation) which
  performs periodic trials to evaluate what codecs work best.  The
  choice of codecs permitted is one of the things tweaked by -1 to -9
  (only odd numbers implemented at the moment).  For example -9 has
  the most number of alternatives for fqzcomp_qual models.

Results
=======

Benchmarks to come...
