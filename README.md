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

ERR174310
---------

Illumina HiSeq 2000 paired end sequencing
Each of the two files compressed independently and results summed together.

This file is also MPEG-G Dataset 01 and has figures reported in their
paper.

| Tool                 | Total(MB)  | Name | Lengths  | Sequence | Qual  |  Elapsed(s) | CPU (s) |
| -------------------- | ---------: | ---: | -------: | -------: | ----: | ----------: | ------: |
| gzip                 |      36652 |      |          |          |       |             |         |
| mpeg-g (no assembly) |      24112 |  376 | (in seq) |    10249 | 13486 |        1193 |    6166 |
| fqzcomp5 -1          |      26652 | 2743 |        0 |    10104 | 13813 |         904 |    1145 |
| fqzcomp5 -3          |      25682 | 1770 |        0 |    10103 | 13808 |        1111 |    1600 |
| fqzcomp5 -5          |      23355 | 1770 |        0 |     9548 | 12037 |        1542 |    5472 |
| fqzcomp5 -7          |      21922 |  668 |        0 |     9288 | 11966 |        2009 |    6576 |
| fqzcomp5 -9          |      21758 |  667 |        0 |     9136 | 11956 |        2435 |    8567 |
| fqzcomp4 n2 s7+b q3  |      21009 |  817 | (in seq) |     8432 | 11759 |     (*)5108 |    813  |


Fqzcomp5's elapsed time is with 4 threads.  GenomSys publication uses 8
threads, with the CPU time being reported as their elapsed 1-thread
benchmark as an estimation.  Both reports used very similar CPUs, but
the I/O systems may be different.  It's clear fqzcomp was I/O bound
for faster methods, and possibly their tool too.

(*) Fqzcomp4 is the smallest, but it lacks random access capability.
Notable is the CPU time which is low, considerably lower than the
elapsed time.  CPU utilisation was bad (about 15%), but this is due to
the bottleneck of reading the input via a zcat pipe.  The elapsed time
should not be considered as meaningful.

The faster rANS/LZP based fqzcomp5 levels (-1 and -3) are still
reasonably competitive and clearly beat the original two fastq.gz
sizes (reported in gzip row).  Once adaptive range coding models are
used the CPU time is considerably slower, but data size drops
considerably.  Level 9 doesn't seem to offer any benefit here over 7.


SRR1238539
----------

IonTorrent variable length WGS data.

This file is also MPEG-G Dataset 11 and has figures reported in their
paper.


| Tool                 | Total(MB)  | Name | Lengths  | Sequence | Qual  |  Elapsed(s) | CPU (s) |
| -------------------- | ---------: | ---: | -------: | -------: | ----: | ----------: | ------: |
| gzip                 |      25000 |      |          |          |       |             |    5608 |
| zstd -1              |      26590 |      |          |          |       |             |    2338 |
| zstd -6              |      24110 |      |          |          |       |             |    3120 |
| zstd -15             |      23205 |      |          |          |       |             |   25038 |
| xz                   |      21324 |      |          |          |       |             |   57910 |
| mpeg-g (no assembly) |      21769 |  205 | (in seq) |     7873 | 13692 |         960 |    5982 |
| fqzcomp5 -1          |      19391 |  818 |      337 |     3742 | 14494 |         378 |     748 |
| fqzcomp5 -3          |      19145 |  583 |      337 |     3733 | 14492 |         441 |     941 |
| fqzcomp5 -5          |      15910 |  583 |      337 |     3733 | 11257 |         713 |    2978 |
| fqzcomp5 -7          |      15353 |  163 |      337 |     3737 | 11116 |         998 |    3484 |
| fqzcomp5 -9          |      15332 |  163 |      337 |     3738 | 11093 |        1063 |    3673 |
| fqzcomp4 n2 s7 b q3  |      20319 |  187 | (in seq) |     6728 | 13241 |        2075 |    3305 |
| colord -balanced     |      17087 |  401 | (in seq) |     4043 | 12643 |        5480 |   23461 |
| colord -ratio        |      16033 |  401 | (in seq) |     2712 | 12921 |       10642 |   37808

As before Fqzcomp's elapsed time is with 4 threads and MPEG-G's is 8 threads.

Fqzcomp's sequence is small here due to the use of LZP and an apparent
genomic ordering in sequences despite being a FASTQ file.  This data
was originally converted by EBI from an SRA file.  Maybe it's possible
SRA has already done some sequence sorting for compression purposes,
but the original read identifiers have been lost by the SRA so it's
not possible to glean the original ordering to verify this.

The gzip size here is as reported by rerunning gzip rather than taking
the input fastq.gz size (it started as .sra anyway), which means we
have times.  Also shown are other generic compression tools.  Xz is
reasonably competitive with light-weight fqz encoding, but is
extremely slow.

Fqzcomp -1 and -3 are fast rANS / LZP / tokenisation strategies.  -5
onwards is using markov modelling and adaptive range coding for
quality values, which markedly reduces the size.  The quality encoding
here is exploiting both previous quality values and neighbouring
sequence bases as context.

Fqzcomp4 performs poorly on this data set, compressing neither the
sequence nor the quality well.

Colord in balanced mode does a reasonable job, but is let down mainly
by the quality compression and the slow speed.  It doesn't have an
IonTorrent profile, but testing on a smaller subset the ONT profile
gave the best ratio so this was used for the full dataset. (On the
first 1 million records, the ONT profile is 10% smaller than the
others for quality values.)  I'm sure it would be more competitive on
sequence compression had the data not apparently already been sorted.
