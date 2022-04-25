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

In the results below I mostly took the standard options for the tools.
It's likely better results are available by changing more parameters,
such as increasing the block sizes.  I did not explore this however
as it's already a large search space.


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
| fqzcomp4 n2 s7+b q3  |      21009 |  817 | (in seq) |     8432 | 11759 |        4901 |    6967 |
| genozip              |      23726 |~1426 |        0 |    ~9771 |~12455 |        5472 |   21620 |
| genozip --best=NO_REF|      22995 |~1426 |        0 |    ~9448 |~12241 |       10376 |   39483 |
| CRAM 3.0             |      25121 | 1100 |  (fixed) |    10139 | 13867 |         908 |    1988 |
| CRAM 3.1 small       |      23231 |  404 |  (fixed) |    10102 | 12719 |         964 |    4740 |
| Spring               |      15368 |  407 |  (fixed) |     2431 | 12523 |       10943 |   35438 |

Fqzcomp5's elapsed time is with 4 threads.  GenomSys publication uses 8
threads, with the CPU time being reported as their elapsed 1-thread
benchmark as an estimation.  Both reports used very similar CPUs, but
the I/O systems may be different.  It's clear fqzcomp was I/O bound
for faster methods, and possibly their tool too.

Fqzcomp4 is the smallest, but it lacks random access capability.

The faster rANS/LZP based fqzcomp5 levels (-1 and -3) are still
reasonably competitive and clearly beat the original two fastq.gz
sizes (reported in gzip row).  Once adaptive range coding models are
used the CPU time is considerably slower, but data size drops
considerably.  Level 9 doesn't seem to offer any benefit here over 7.

Genozip-13.0.5 was ran in the default mode (as "best" needs a
reference).  Sizes are approximate due to coarse rounding in the
genocat --stats output.

CRAM sizes are a little different to fqzcomp here as the command used was

```
   samtools import -N -1 ERR174310_1.fastq.gz -2 ERR174310_2.fastq.gz
```

This interleaves the two files rather than compressing each separately
and aggregating results, which in turn halves the size of the read
name encoding due to deduplication.  It also performs reasonably well
with CRAM 3.1, but has poorer quality compression.  This is perhaps
due to finer-grained random access (unverified).   CRAM 3.1 lacks an
adaptive sequence model, so is behind fqzcomp on sequence
compression.

Spring also interleaves both files.  Spring is doing a local read
clustering and reordering process to help compress the sequence data.
This is both CPU and memory intensive (about 21GB on this data set).
While it can achieve good results, I would argue that the resources
would be best spent on sequence alignment and/or a proper denovo
assembly (or maybe a mixture of both, with assembly on the data that
doesn't map to identify the large insertions, organelles and
contaminants).  This is probably a significant step up again in
resources, but the end result is far more useful to the user.

I view FASTQ largely as an interim data format - somewhere between the
raw instrument data and the end product (assembled BAM/CRAM and/or
VCFs).  It doesn't make a great deal of sense to me to spend a huge
amount of CPU on FASTQ compression alone.


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
| colord -ratio        |      16033 |  401 | (in seq) |     2712 | 12921 |       10642 |   37808 |
| genozip              |      18926 | ~480 | (in seq) |    ~5261 |~13207 |        3776 |   14682 |
| genozip --best=NO_REF|      18511 | ~460 | (in seq) |    ~5154 |~12885 |        6509 |   25288 |
| CRAM 3.0             |      18340 |  316 |      209 |     3320 | 14490 |         370 |    1719 |
| CRAM 3.1 small       |      17302 |   ?  |      183 |     3340 | 13774 |        1239 |    5169 |

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

Genozip-13.0.5 was ran in the default mode and "--best=NO_REF".

CRAM on this data set performs well at lighter levels in CRAM 3.0,
again due to the unusual nature of the data already being
clustered by sequence identity.  CRAM 3.1 currently lacks the ability
to use sequence bases in the FQZ-qual model, so loses out considerably
at higher compression levels.


ERR2442595
----------

This is a relatively small RNASeq dataset sequenced on Illumina
NovaSeq 6000.  The quality on this instrument is binned to just 4
discrete values.  The nature of the RNAseq experiment gives an
effective small portion of genome covered, which aids sequence
compression.

| Tool                 | Total(MB)  | Name | Lengths  | Sequence | Qual  |  Elapsed(s) | CPU (s) |
| -------------------- | ---------: | ---: | -------: | -------: | ----: | ----------: | ------: |
| gzip                 |       3852 |      |          |          |       |             |         |
| fqzcomp5 -1          |       2625 |  132 |        0 |     2076 |   417 |         118 |     151 |
| fqzcomp5 -3          |       2538 |   60 |        0 |     2067 |   411 |         121 |     204 |
| fqzcomp5 -5          |       1817 |   60 |        0 |     1361 |   396 |         236 |     873 |
| fqzcomp5 -7          |       1406 |    0 |        0 |     1014 |   392 |         407 |    1359 |
| fqzcomp5 -9          |       1193 |    0 |        0 |      801 |   392 |         696 |    2316 |
| fqzcomp4 n2 s7 b q3  |       1038 |    0 | (in seq) |      636 |   402 |         973 |    1022 |
| spring               |        645 |    0 |        0 |      233 |   410 |         682 |    1970 |
| genozip --best=NO_REF|       1464 |   60 |        0 |      988 |   416 |         702 |    1926 |


Spring is included here to show the impact of sequence reordering and
clustering.  It significantly reduces compressed sequence size over
the fqzcomp statistical model.  On this data set it's doesn't have a
large memory and CPU usage, although this can become an issue with
bigger datasets as seen above

