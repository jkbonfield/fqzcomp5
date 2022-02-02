/* Tests for fqz codec */
/*
 * Copyright (c) 2019,2020,2022 Genome Research Ltd.
 * Author(s): James Bonfield
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *    1. Redistributions of source code must retain the above copyright notice,
 *       this list of conditions and the following disclaimer.
 *
 *    2. Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *    3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
 *       Institute nor the names of its contributors may be used to endorse
 *       or promote products derived from this software without specific
 *       prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH
 * LTD OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <fcntl.h>
#include <ctype.h>
#include <limits.h>

#include "htscodecs/fqzcomp_qual.h"
#include "htscodecs/varint.h"

#ifndef MAX_REC
#define MAX_REC 1000000
#endif

#ifndef MAX_SEQ
#  define MAX_SEQ 100000
#endif

#ifndef MIN
#  define MIN(a,b) ((a)<(b)?(a):(b))
#  define MAX(a,b) ((a)>(b)?(a):(b))
#endif

static fqz_slice fixed_slice = {0};

fqz_slice *fake_slice(size_t buf_len, int *len, int *r2, int *sel, int nlen,
		      unsigned char *seq, int *seq_idx) {
    fixed_slice.num_records = (nlen == 1 && len)
	? (buf_len+len[0]-1) / len[0] : nlen;
    assert(fixed_slice.num_records <= MAX_REC);
    int i;
    if (!fixed_slice.len)
	fixed_slice.len = malloc(MAX_REC * sizeof(*fixed_slice.len));
    if (!fixed_slice.flags)
	fixed_slice.flags = malloc(MAX_REC * sizeof(*fixed_slice.flags));
    if (!fixed_slice.seq)
	fixed_slice.seq = malloc(MAX_REC * sizeof(*fixed_slice.seq));
    for (i = 0; i < fixed_slice.num_records; i++) {
	int idx = i < nlen ? i : nlen-1;
	fixed_slice.len[i] = len ? len[idx] : 0;
	fixed_slice.flags[i] = r2 ? r2[idx]*FQZ_FREAD2 : 0;
	fixed_slice.flags[i] |= sel ? (sel[idx]<<16) : 0;
	fixed_slice.seq[i] = seq_idx ? seq + seq_idx[i] : NULL;
    }

    return &fixed_slice;
}

typedef struct {
    int num_records;
    char *buf;
    char **name; // pointers into buf
    char **seq;  // pointers into buf
    char **qual; // pointers into buf
    int *len;    // sequence length
    int *flag;   // READ1/READ2 parsed from name
} fastq;

// Load a line into &fq->buf[i], alloced to abuf.
size_t get_line(FILE *fp, fastq *fq, size_t i, size_t *abuf) {
    for(;;) {
	if (!fgets(&fq->buf[i], *abuf-i, fp))
	    return 0;
	size_t l = strlen(&fq->buf[i]);
	i += l;
	if (fq->buf[i-1] == '\n')
	    break;

	if (i >= *abuf) {
	    *abuf = *abuf*1.5 + 10000;
	    fq->buf = realloc(fq->buf, *abuf);
	    if (!fq->buf)
		return 0;
	}
    }

    fq->buf[i-1] = 0;
    return i;
}

fastq *fastq_alloc(int nr) {
    fastq *fq = calloc(1, sizeof(*fq));

    fq->num_records = nr;
    fq->name = calloc(nr, sizeof(char *));
    fq->seq  = calloc(nr, sizeof(char *));
    fq->qual = calloc(nr, sizeof(char *));
    fq->len  = calloc(nr, sizeof(int));
    fq->flag = calloc(nr, sizeof(int));

    return fq;
}

void fastq_free(fastq *fq) {
    if (!fq)
	return;
    free(fq->buf);
    free(fq->name);
    free(fq->seq);
    free(fq->qual);
    free(fq->len);
    free(fq->flag);
    free(fq);
}

// Debug function
void fastq_dump(fastq *fq) {
    int i;
    for (i = 0; i < fq->num_records; i++) {
	printf("@%s\n%.*s\n+\n%.*s\n",
	       fq->name[i], fq->len[i], fq->seq[i], fq->len[i], fq->qual[i]);
    }
}

#define goto if (fprintf(stderr, "ERR %s:%d\n", __FILE__, __LINE__)) goto
fastq *load_seqs(FILE *in, int blk_size) {
    fastq *fq = calloc(1, sizeof(*fq));
    if (!fq)
	goto err;
    size_t abuf = blk_size+10000;
    char *buf = fq->buf = malloc(abuf);
    if (!buf)
	goto err;
    int nseq = 0, aseq = 0;

    // Load initial tranch of data
    size_t n = fread(buf, 1, blk_size, in);
    if (n <= 0) {
	if (feof(in))
	    return fq;
	else
	    goto err;
    }

    // Parse lines we've loaded so far, in sets of 4 (simple FASTQ only)
    size_t i = 0;
    char *name, *seq, *qual, *last_name = NULL;
    int line = 0;

    while (i < n) {
	int j = i;

	// @name
	if (j < n && buf[j] != '@')
	    goto err;
	name = &buf[++j];
	while (j < n && buf[j] != '\n')
	    j++;
	if (j >= n || buf[j] != '\n')
	    break;
	buf[j] = 0;
	line++;

	// seq
	seq = &buf[++j];
	while (j < n && buf[j] != '\n')
	    j++;
	if (j >= n || buf[j] != '\n')
	    break;
	buf[j] = 0;
	line++;

	// +(name)
	if (++j < n && buf[j] != '+')
	    goto err;
	while (j < n && buf[j] != '\n')
	    j++;
	if (j >= n || buf[j] != '\n')
	    break;
	line++;

	// qual
	qual = &buf[++j];
	while (j < n && buf[j] != '\n')
	    j++;
	if (j >= n || buf[j] != '\n')
	    break;
	buf[j] = 0;

	int len;
	if ((len = strlen(seq)) != strlen(qual))
	    goto err;
	line++;

	if (nseq+1 >= aseq) {
	    aseq = aseq*1.5 + 1000;
	    fq->name = realloc(fq->name, aseq * sizeof(char *));
	    fq->seq  = realloc(fq->seq,  aseq * sizeof(char *));
	    fq->qual = realloc(fq->qual, aseq * sizeof(char *));
	    fq->len  = realloc(fq->len,  aseq * sizeof(int));
	    fq->flag = realloc(fq->flag, aseq * sizeof(int));
	    if (!fq->name || !fq->seq || !fq->qual || !fq->len || !fq->flag)
		goto err;
	}
	fq->name[nseq] = name;
	fq->seq [nseq] = seq;
	fq->qual[nseq] = qual;
	fq->len [nseq] = len;

	// name/1 and name/2 for identifying read-pairs,
	// or n1 n1 n2 n2 n3 n3 pairing detection.
	len = strlen(name);
	int flag = 0;
	if (len > 2 && name[len-1] == '2' && name[len-2] == '/')
	    flag = FQZ_FREAD2;
	if (last_name && strcmp(name, last_name))
	    flag = FQZ_FREAD2;
	fq->flag[nseq] = flag;
	last_name = name;
	nseq++;

	i = ++j;
	if (i < n)
	    line = 0;
    }
    fq->num_records = nseq;

    // Load remaining partial fastq record
    i = n;
    switch (line) {
    case 0: // @name
	if (!(i = get_line(in, fq, i, &abuf)))
	    goto err;
	seq = &fq->buf[i];
	// fall through
	
    case 1: // seq
	if (!(i = get_line(in, fq, i, &abuf)))
	    goto err;
	// fall through

    case 2: // +(name)
	if (!(i = get_line(in, fq, i, &abuf)))
	    goto err;
	qual = &fq->buf[i];
	// fall through

    case 3: // qual
	if (!(i = get_line(in, fq, i, &abuf)))
	    goto err;
	// fall through

	int len;
	if ((len = strlen(seq)) != strlen(qual))
	    goto err;

	fq->name[nseq] = name;
	fq->seq [nseq] = seq;
	fq->qual[nseq] = qual;
	fq->len [nseq] = len;

	len = strlen(name);
	int flag = 0;
	if (len > 2 && name[len-1] == '2' && name[len-2] == '/')
	    flag = FQZ_FREAD2;
	if (last_name && strcmp(name, last_name))
	    flag = FQZ_FREAD2;
	fq->flag[nseq] = flag;
	last_name = name;


	fq->num_records = ++nseq;

    case 4: // complete
	break;
    }

    return fq;

 err:
    fprintf(stderr, "Failed to load fastq input\n");
    fastq_free(fq);

    return NULL;
}

static uint64_t manual_strats[10] = {0};
static int manual_nstrat = 0;

/*
 * Manually specified strategies held in global manual_strats[].
 */
static inline
int fqz_manual_parameters(fqz_gparams *gp,
			  fqz_slice *s,
			  unsigned char *in,
			  size_t in_size) {
    int i, p;
    int dsqr[] = {
	0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5,
	5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
	6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
    };

    gp->vers = FQZ_VERS;
    gp->nparam = manual_nstrat;
    gp->gflags = GFLAG_MULTI_PARAM | GFLAG_HAVE_STAB;
    for (i = 0; i < 256; i++)
	gp->stab[i] = 0;

    // Fill these out later
    gp->max_sel = 0;
    gp->max_sym = 0;
    gp->p = malloc(gp->nparam * sizeof(*gp->p));

    for (p = 0; p < gp->nparam; p++) {
	fqz_param *pm = &gp->p[p];
	uint64_t st = manual_strats[p];

	pm->boff   = st & 15; st >>= 4;
	pm->bloc   = st & 15; st >>= 4;
	pm->bbits  = st & 15; st >>= 4;
	pm->do_qa  = st & 15; st >>= 4;
	pm->do_r2  = st & 15; st >>= 4;
	pm->dloc   = st & 15; st >>= 4;
	pm->ploc   = st & 15; st >>= 4;
	pm->sloc   = st & 15; st >>= 4;
	pm->qloc   = st & 15; st >>= 4;
	pm->dshift = st & 15; st >>= 4;
	pm->dbits  = st & 15; st >>= 4;
	pm->pshift = st & 15; st >>= 4;
	pm->pbits  = st & 15; st >>= 4;
	pm->qshift = st & 15; st >>= 4;
	pm->qbits  = st & 15; st >>= 4;

	// Gather some stats, as per qual_stats func.
	// r in rec count.
	// i = index to in[]
	// j = index within this rec
	uint32_t qhist[256] = {0};

	// qual stats for seqs using this parameter only
	fqz_qual_stats(s, in, in_size, pm, qhist, p);
	int max_sel = pm->max_sel;

	// Update max_sel running total. Eg with 4 sub-params:
	//
	// sel    param no.   => new
	// 0      0              0
	// 0/1    1              1,2
	// 0/1    2              3,4
	// 0      3              5
	for (i = gp->max_sel; i < gp->max_sel + max_sel+1; i++)
	    gp->stab[i] = p;
	gp->max_sel += max_sel+1;

	pm->fixed_len = pm->fixed_len > 0;
	pm->use_qtab = 0;  // unused by current encoder
	pm->store_qmap = pm->nsym <= 8;

	// Adjust parameters based on quality stats.
	// FIXME:  dup from fqz_pick_parameters.
	for (i = 0; i < sizeof(dsqr)/sizeof(*dsqr); i++)
	    if (dsqr[i] > (1<<pm->dbits)-1)
		dsqr[i] = (1<<pm->dbits)-1;

	if (pm->store_qmap) {
	    int j;
	    for (i = j = 0; i < 256; i++)
		if (qhist[i])
		    pm->qmap[i] = j++;
		else
		    pm->qmap[i] = INT_MAX;
	    pm->max_sym = pm->nsym;
	} else {
	    pm->nsym = 255;
	    for (i = 0; i < 256; i++)
		pm->qmap[i] = i;
	}
	if (gp->max_sym < pm->max_sym)
	    gp->max_sym = pm->max_sym;

	// Produce ptab from pshift.
	if (pm->qbits) {
	    for (i = 0; i < 256; i++) {
		pm->qtab[i] = i; // 1:1
		//pm->qtab[i] = (1<<pm->qshift)-i; // 1:1
		//pm->qtab[i] = i/4;
		//pm->qtab[i] = ((1<<pm->qshift)-i)/3;

		// Alternative mappings:
		//qtab[i] = i > 30 ? MIN(max_sym,i)-15 : i/2;  // eg for 9827 BAM
	    }
#if 0
	    // qtab for PacBio CCS data; saves 3%
	    for (i='~'-33; i<256; i++) {
		pm->qtab[i] = 24+i-('~'-33);
	    }

	    int x = 0;
	    for (i = 0; i < 1; i++)
		pm->qtab[i] = x,x++;
	    for (;i < '~'-33; i++)
		//pm->qtab[i] = x,x+=(i%4==0);
		pm->qtab[i] = x,x+=(i%4==0);
	    x++;
	    for (;i <= '~'-33; i++)
		pm->qtab[i] = x,x++;
	    for (;i < 256; i++)
		pm->qtab[i] = x,x+=(i%4==0);
#endif

//	    for (i = 0; i < 128; i++) {
//		for (x = i+1; x < 128; x++) {
//		    if (pm->qtab[i] != pm->qtab[x])
//			break;
//		}
//		x--;
//		if (i==x)
//		    fprintf(stderr, "%d:%d ", pm->qtab[i], i);
//		else {
//		    fprintf(stderr, "%d:%d-%d ", pm->qtab[i], i, x);
//		    i=x;
//		}
//	    }
//	    fprintf(stderr, "\n");

	    //pm->qtab['~'-33]=32;

	    // pm->use_qtab = 1;
//	    for (i = 0; i <= 2 ; i++) pm->qtab[i] = 0;
//	    for (     ; i <= 12; i++) pm->qtab[i] = 1;
//	    for (     ; i <= 18; i++) pm->qtab[i] = 2;
//	    for (     ; i <= 36; i++) pm->qtab[i] = 3;
	}
	//pm->use_qtab = 0;
	pm->qmask = (1<<pm->qbits)-1;

	if (pm->pbits) {
	    for (i = 0; i < 1024; i++)
		pm->ptab[i] = MIN((1<<pm->pbits)-1, i>>pm->pshift);

//	    for (i = 0; i < 1024; i++)
//		pm->ptab[i] = MIN((1<<pm->pbits)-1, i < 10 ? i : 10 + i/3);

	    // Alternatively via analysis of quality distributions we
	    // may select a bunch of positions that are special and
	    // have a non-uniform ptab[].
	    // Manual experimentation on a NovaSeq run saved 2.8% here.
	}

	if (pm->dbits) {
	    for (i = 0; i < 256; i++)
		pm->dtab[i] = dsqr[MIN(sizeof(dsqr)/sizeof(*dsqr)-1, i>>pm->dshift)];
	}

	pm->use_ptab = (pm->pbits > 0);
	pm->use_dtab = (pm->dbits > 0);

	pm->pflags =
	    (pm->use_qtab   ?PFLAG_HAVE_QTAB :0)|
	    (pm->use_dtab   ?PFLAG_HAVE_DTAB :0)|
	    (pm->use_ptab   ?PFLAG_HAVE_PTAB :0)|
	    (pm->do_sel     ?PFLAG_DO_SEL    :0)|
	    (pm->fixed_len  ?PFLAG_DO_LEN    :0)|
	    (pm->do_dedup   ?PFLAG_DO_DEDUP  :0)|
	    (pm->store_qmap ?PFLAG_HAVE_QMAP :0);
    }

    for (i = gp->max_sel; i < 256; i++)
	gp->stab[i] = gp->stab[gp->max_sel-1];

    return 0;
}

#define BS 1024*1024
static unsigned char *load(char *fn, size_t *lenp) {
    unsigned char *data = NULL;
    uint64_t dsize = 0;
    uint64_t dcurr = 0;
    signed int len;

    //build_rcp_freq();

#ifndef _O_BINARY
#define _O_BINARY 0
#endif

    int fd = open(fn, O_RDONLY | _O_BINARY);
    if (fd < 0) {
	perror(fn);
	return NULL;
    }

    do {
	if (dsize - dcurr < BS) {
	    dsize = dsize ? dsize * 2 : BS;
	    data = realloc(data, dsize);
	}

	len = read(fd, data + dcurr, BS);
	if (len > 0)
	    dcurr += len;
    } while (len > 0);

    if (len == -1) {
	perror("read");
    }
    close(fd);

    *lenp = dcurr;
    return data;
}

#define BLK_SIZE 300*1000000
//#define BLK_SIZE 10*1000000

int count_lines(unsigned char *in, size_t len) {
    size_t i;
    int lines = 0;

    for (i = 0; i < len; i++)
	if (in[i] == '\n')
	    lines++;

    return lines;
}

// QUAL [is_read2 [selector]]
void parse_lines(unsigned char *in, size_t len,
		 int *rec_len, int *rec_r2, int *rec_sel,
		 size_t *new_len) {
    size_t i, j, start;
    int rec = 0;

    for (start = i = j = 0; i < len; i++) {
	if (in[i] == '\n' || in[i] == ' ' || in[i] == '\t') {
	    rec_len[rec] = i-start;

	    // Read2 marker
	    while (i < len && in[i] != '\n' && isspace(in[i]))
		i++;

	    if (in[i] != '\n')
		rec_r2[rec] = atoi((char *)&in[i]);
	    else
		rec_r2[rec] = 0;

	    while (i < len && !isspace(in[i]))
		i++;

	    // selector
	    while (i < len && in[i] != '\n' && isspace(in[i]))
		i++;

	    if (in[i] != '\n')
		rec_sel[rec] = atoi((char *)&in[i]);
	    else
		rec_sel[rec] = 0;

	    while (i < len && in[i] != '\n')
		i++;

	    start = i+1;
	    rec++;
	} else {
	    in[j++] = in[i]-33; // ASCII phred to qual
	}
    }

    *new_len = j;
}

void parse_seq(unsigned char *seq, int *seq_idx, int nrec) {
    int i, j;
    if (!seq)
	return;

    for (i = j = 0; i < nrec; i++) {
	seq_idx[i] = j;
	while (seq[j] != '\n')
	    j++;
	j++;
    }
}


int main(int argc, char **argv) {
    unsigned char *in, *out, *seq;
    size_t in_len, out_len, seq_len;
    int decomp = 0, vers = 4;  // CRAM version 4.0 (4) or 3.1 (3)
    int strat = 0, raw = 0;
    fqz_gparams *gp = NULL, gp_local;
    int blk_size = BLK_SIZE; // MAX
    FILE *in_fp, *out_fp;

#ifdef _WIN32
        _setmode(_fileno(stdin),  _O_BINARY);
        _setmode(_fileno(stdout), _O_BINARY);
#endif

    extern char *optarg;
    extern int optind;
    int opt;

    while ((opt = getopt(argc, argv, "ds:s:b:rx:")) != -1) {
	switch (opt) {
	case 'd':
	    decomp = 1;
	    break;

	case 'b':
	    blk_size = atoi(optarg);
	    if (blk_size > BLK_SIZE)
		blk_size = BLK_SIZE;
	    break;

	case 's':
	    strat = atoi(optarg);
	    break;

	case 'x': {
	    // Hex digits are:
	    // qbits  qshift
	    // pbits  pshift
	    // dbits  dshift
	    // qloc   sloc
	    // ploc   dloc
	    // do_r2  do_qavg
	    //
	    // Examples: -x 0x5570000d6e14 q40+dir =  3473340
	    //           -x 0x8252120e8d04 q4      =  724989
	    uint64_t x = strtol(optarg, NULL, 0);
	    manual_strats[manual_nstrat++] = x;

	    gp = &gp_local;
	    break;
	}
	case 'r':
	    raw = 1;
	    break;
	}
    }

    in_fp = optind < argc ? fopen(argv[optind], "r") : stdin;
    out_fp = ++optind < argc ? fopen(argv[optind], "wb") : stdout;

    // Block based, for arbitrary sizes of input
    if (decomp) {
	for (;;) {
	    // Load next compressed block
	    int nb, nr, c_len, i;
	    if (fread(&nb, 1, 4, in_fp) != 4)
		break;
	    if (fread(&nr, 1, 4, in_fp) != 4)
		break;
	    if (fread(&c_len, 1, 4, in_fp) != 4)
		break;
	    fastq *fq = fastq_alloc(nr);

	    // ----------
	    // Qual
	    char *comp = malloc(c_len);
	    if (fread(comp, 1, c_len, in_fp) != c_len)
		break;

	    int *lengths = malloc(nr * sizeof(lengths));
	    out = fqz_decompress((char *)comp, c_len, &out_len,
				 lengths, nr, NULL);

	    // ----------
	    // Convert back to fastq
	    char *cp = out;
	    for (i = 0; i < nr; i++) {
		printf("%.*s\n", lengths[i], cp);
		cp += lengths[i];
	    }

	    free(lengths);
	    free(out);
	    free(fq);
	}

//	fastq *fq = NULL;
//
//	if (seq) {
//	    int nlines = count_lines(seq, seq_len);
//	    int *seq_idx = calloc(nlines, sizeof(*seq_idx));
//	    parse_seq(seq, seq_idx, nlines);
//
//	    //s = fake_slice(0, NULL, NULL, NULL, nlines, seq, seq_idx);
//	}
//
//	unsigned char *in2 = in;
//	while (in_len > 0) {
//	    // Read sizes as 32-bit
//	    size_t in2_len, out_len;
//	    if (raw) {
//		uint32_t u32;
//		var_get_u32(in2, in2+in_len, &u32);
//		out_len = u32;
//		in2_len = in_len;
//	    } else {
//		out_len = *(uint32_t *)in2;  in2 += 4;
//		in2_len = *(uint32_t *)in2;  in2 += 4;
//	    }
//
//	    fprintf(stderr, "out_len %ld, in_len %ld\n", (long)out_len, (long)in2_len);
//
//	    int *lengths = malloc(MAX_REC * sizeof(int));
//	    out = (unsigned char *)fqz_decompress((char *)in2, in_len-(raw?0:8), &out_len, lengths, MAX_REC, s);
//	    if (!out) {
//		fprintf(stderr, "Failed to decompress\n");
//		return 1;
//	    }
//
//	    // Convert from binary back to ASCII with newlines
//	    int i = 0, j = 0;
//	    while (j < out_len) {
//		int k;
//		char seq[MAX_SEQ];
//		for (k = 0; k < lengths[i]; k++)
//		    seq[k] = out[j+k]+33;
//		seq[k] = 0;
//		puts(seq);
//		j += lengths[i++];
//	    }
//	    free(out);
//	    in2 += in2_len;
//	    in_len -= in2_len+(raw?0:8);
//
//	    free(lengths);
//
//	    break; // One cycle only until we fix blocking to be \n based
//	}
    } else {
	// Load the next batch of fastq entries
	for(;;) {
	    fastq *fq = load_seqs(in_fp, blk_size);
	    if (!fq)
		exit(1);

	    if (!fq->num_records) {
		fastq_free(fq);
		break;
	    }

	    //fastq_dump(fq);

	    //----------
	    // Names: tok3

	    //----------
	    // Seq: rans or statistical modelling

	    //----------
	    // Qual: rans or fqz
	    // Convert fastq struct to fqz_slice for context
	    fqz_slice *s = malloc(fq->num_records * sizeof(*s));
	    s->num_records = fq->num_records;
	    s->len = fq->len;
	    s->flags = fq->flag;
	    s->seq = (unsigned char **)fq->seq;

	    // Concatenate qualities together into a single block.
	    int nb = 0, i, j;
	    for (i = 0; i < fq->num_records; i++)
		nb += fq->len[i];
	    uint8_t *qual = malloc(nb);
	    for (i = j = 0; i < fq->num_records; j += fq->len[i++])
		memcpy(&qual[j], fq->qual[i], fq->len[i]);
	    assert(nb == j);

	    out = fqz_compress(vers, s, (char *)qual, nb, &out_len, strat, gp);
	    fwrite(&nb, 1, 4, out_fp);      // FIXME: use var_put_u32
	    fwrite(&fq->num_records, 1, 4, out_fp);
	    fwrite(&out_len, 1, 4, out_fp); // FIXME: use var_put_u32
	    fwrite(out, 1, out_len, out_fp);
	    free(out);

	    fastq_free(fq);
	    free(s);
	    free(qual);
	}
    }

    fclose(in_fp);
    fclose(out_fp);

    return 0;
}
