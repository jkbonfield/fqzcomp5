// TODO
// - Tokenise name with different alphabets.
//   See ~/scratch/data/enano_vir/ERR2708432.fastq

// - Also permit old name encoding (fqzcomp -n1).
//   Better on above data. Why?
//   (fqz -n1: 912989, -n2: 1340223, paq8: 891928

// - Separate name vs comment field encoding.

// - Fqzcomp_qual for sequences
//   Or optional bi-stranded analysis

// - Removal of lengths from fqzqual stream

// - Entropy encoding of read length stream

// - Exploration of multiple rans options (o4, o5, o133, o197?)

// - Reuse of memory buffers for sped

// - Threading!

// - Increase sizes of fqzcomp contexts?  Not so useful for CRAM, but maybe
//   it's still beneficial to go beyond 16-bit.

// - Improve fqzcomp tables to permit any mapping rather than monotonic.

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

#include "htscodecs/varint.h"
#include "htscodecs/fqzcomp_qual.h"
#include "htscodecs/tokenise_name3.h"
#include "htscodecs/rANS_static4x16.h"
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

typedef struct {
    int num_records;                 // number of fastq entries
    char *name_buf;                  // concatenated names, \0 separator
    char *seq_buf;                   // concatenated seq,   no separator
    char *qual_buf;                  // concatenated qual,  no separator
    int *name;                       // index into name_buf
    int *seq;                        // index into seq_buf
    int *qual;                       // index into qual_buf
    int *len;                        // sequence length
    int *flag;                       // READ1/READ2 parsed from name
    int name_len, seq_len, qual_len; // used size of _buf above
    int name_sz,  seq_sz,  qual_sz;  // alloced size of _buf above
    int fixed_len;                   // length of each seq, 0 if not fixed
} fastq;

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
    free(fq->name_buf);
    free(fq->seq_buf);
    free(fq->qual_buf);
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
	       fq->name_buf + fq->name[i],
	       fq->len[i], fq->seq_buf + fq->seq[i],
	       fq->len[i], fq->qual_buf + fq->qual[i]);
    }
}

#define goto if (fprintf(stderr, "ERR %s:%d\n", __FILE__, __LINE__)) goto
fastq *load_seqs(FILE *in, int blk_size) {
    fastq *fq = calloc(1, sizeof(*fq));
    if (!fq)
	goto err;
    size_t name_sz = blk_size/4;
    size_t seq_sz  = blk_size/2;
    size_t qual_sz = blk_size/2;
    char *name_buf = fq->name_buf = malloc(name_sz);
    char *seq_buf  = fq->seq_buf  = malloc(seq_sz);
    char *qual_buf = fq->qual_buf = malloc(qual_sz);
    if (!name_buf || !seq_buf || !qual_buf)
	goto err;
    int c, i = 0, nr = 0, ar = 0;
    int last_name = -1;
    fq->fixed_len = -1;

    int name_i = 0, seq_i = 0, qual_i = 0;
    while (i < blk_size) {
	if (nr >= ar) {
	    ar = ar*1.5 + 10000;
	    fq->name = realloc(fq->name, ar*sizeof(char *));
	    fq->seq  = realloc(fq->seq , ar*sizeof(char *));
	    fq->qual = realloc(fq->qual, ar*sizeof(char *));
	    fq->len  = realloc(fq->len,  ar*sizeof(int));
	    fq->flag = realloc(fq->flag, ar*sizeof(int));
	}

	// @name
	fq->name[nr] = name_i;
	if ((c = getc(in)) != '@') {
	    if (c == EOF)
		break;
	    else
		goto err;
	}
	while ((c = getc(in)) != EOF && c != '\n') {
	    if (name_i >= name_sz) {
		name_sz = name_sz * 1.5 + 1000;
		name_buf = fq->name_buf = realloc(fq->name_buf, name_sz);
	    }
	    name_buf[name_i++] = c;
	}
	if (c == EOF)
	    goto err;
	name_buf[name_i++] = 0;
	i += name_i - fq->name[nr];

	int flag = 0;
	if (name_i > 2 &&
	    name_buf[name_i-1] == '2' &&
	    name_buf[name_i-2] == '/')
	    flag = FQZ_FREAD2;
	if (last_name >= 0 &&
	    strcmp(fq->name_buf + fq->name[nr], fq->name_buf + last_name))
	    flag = FQZ_FREAD2;
	fq->flag[nr] = flag;
	last_name = fq->name[nr];

	// seq
	fq->seq[nr] = seq_i;
	int len = seq_i;
	while ((c = getc(in)) != EOF && c != '\n') {
	    if (seq_i >= seq_sz) {
		// very unlikely given blk_size/2 starting point,
		// but not impossible.
		seq_sz = seq_sz * 1.5 + 1000;
		seq_buf = fq->seq_buf = realloc(fq->seq_buf, seq_sz);
	    }
	    seq_buf[seq_i++] = c;
	}
	if (c == EOF)
	    goto err;
	fq->len[nr] = seq_i - len;
	//seq_buf[seq_i++] = '\n'; // use instead of length terminator?
	i += seq_i - fq->seq[nr];

	if (fq->fixed_len == -1)
	    fq->fixed_len = fq->len[nr];
	else if (fq->fixed_len > 0)
	    if (fq->fixed_len != fq->len[nr])
		fq->fixed_len = 0;

	// +(name)
	if ((c = getc(in)) != '+')
	    goto err;
	while ((c = getc(in)) != EOF && c != '\n')
	    i++;
	if (c == EOF)
	    goto err;
	i++;

	// qual
	fq->qual[nr] = qual_i;
	len = qual_i;
	while ((c = getc(in)) != EOF && c != '\n') {
	    if (qual_i >= qual_sz) {
		// very unlikely given blk_size/2 starting point
		qual_sz = qual_sz * 1.5 + 1000;
		qual_buf = fq->qual_buf = realloc(fq->qual_buf, qual_sz);
	    }
	    qual_buf[qual_i++] = c-33;
	}
	if (c == EOF)
	    goto err;
	if (fq->len[nr] != qual_i - len)
	    goto err;
	i += qual_i - fq->qual[nr];

	nr++;
    }
    fq->name_len = name_i;
    fq->seq_len  = seq_i;
    fq->qual_len = qual_i;
    fq->num_records = nr;

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

#define NSYM 256
#define STEP 8
// FIXME: a dedicated base encoder would be much faster here
#include "htscodecs/c_simple_model.h"
#undef NSYM
#undef STEP

// Why is NSYM 5 so much poorer than NSYM 4?
// Do all those remainder probabilities really add up to be that significant?
#define NSYM 4
#define STEP 1
#include "htscodecs/c_simple_model.h"
// An order-N arithmetic encoder, dedicated to sequence contexts.
char *encode_seq(unsigned char *in,  unsigned int in_size,
		 int ctx_size,
		 /*unsigned char *out,*/ unsigned int *out_size) {
    char *out = malloc(in_size + 100);
    if (!out)
	return NULL;

    const int msize = 1<<(2*ctx_size);
    const int mask = msize-1;

    SIMPLE_MODEL(NSYM,_) *seq_model = malloc(msize * sizeof(*seq_model));
    int i;

    // Do histogram to get observed values.
    // Then set m to max number of elements in histogram.
    int m = NSYM;

    for (i = 0; i < msize; i++)
	SIMPLE_MODEL(NSYM,_init)(&seq_model[i], m);

    SIMPLE_MODEL(256,_) run_len4;
    SIMPLE_MODEL(256,_init)(&run_len4, 256);
    SIMPLE_MODEL(256,_) run_lenN;
    SIMPLE_MODEL(256,_init)(&run_lenN, 256);
    SIMPLE_MODEL(256,_) literal;
    SIMPLE_MODEL(256,_init)(&literal, 256);

    RangeCoder rc;
    RC_SetOutput(&rc, (char *)out);
    RC_StartEncode(&rc);

    int last, last2;
    /* Corresponds to a 12-mer word that doesn't occur in human genome. */
    last  = 0x007616c7 & mask;
    //last2 = (0x2c6b62ff >> (32 - 2*ctx)) & mask; // both strands mode

    int L[256];
    for (int i = 0; i < 256; i++)
	L[i] = 4; // N
    L['A'] = L['a'] = 0;
    L['C'] = L['c'] = 1;
    L['G'] = L['g'] = 2;
    L['T'] = L['t'] = 3;

    // TODO: add an escape mechanism for within ACGT.
    // Like RLE it could be a count of ACGTs, followed by a
    // count of non-ACGTs and a different simpler model.
#if NSYM>4
    for (int i = 0; i < in_size; i++) {
	unsigned char b = L[in[i]];
	SIMPLE_MODEL(NSYM, _encodeSymbol)(&seq_model[last], &rc, b);
	last = ((last<<2) + b) & mask;
    }
#else
    for (int i = 0; i < in_size; ) {//i++) {
	// ACGT symbols
	int j;
	for (j = i; j < in_size; j++) {
	    if (L[in[j]] == 4)
		break;
	}
	int run = j-i, r2 = run;
	for (;;) {
	    SIMPLE_MODEL(256, _encodeSymbol)(&run_len4, &rc, MIN(255, r2));
	    //fprintf(stderr, "run4 %d\n", MIN(255, r2));
	    if (r2 >= 255)
		r2 -= 255;
	    else
		break;
	}

	for (j = 0; j < run; j++) {
	    unsigned char b = L[in[i+j]];
	    SIMPLE_MODEL(NSYM, _encodeSymbol)(&seq_model[last], &rc, b);
	    last = ((last<<2) + b) & mask;

//	    if (1) { // both strands
//		int b2 = last2 & 3;
//		last2 = last2/4 + ((3-b) << (2*NS-2));
//		_mm_prefetch((const char *)&model_seq8[last2], _MM_HINT_T0);
//		SIMPLE_MODEL(NSYM, _updateSymbol)(b2);
//	    }
	}
	i += run;
	if (i >= in_size)
	    break;

	// non-ACGT symbols
	for (j = i; j < in_size; j++) {
	    if (L[in[j]] != 4)
		break;
	}
	run = j-i; r2 = run;
	for (;;) {
	    SIMPLE_MODEL(256, _encodeSymbol)(&run_lenN, &rc, MIN(255, r2));
	    //fprintf(stderr, "runN %d\n", MIN(255, r2));
	    if (r2 >= 255)
		r2 -= 255;
	    else
		break;
	}

	for (j = 0; j < run; j++)
	    SIMPLE_MODEL(256, _encodeSymbol)(&literal, &rc, in[i+j]);
	i += run;
    }
#endif

    RC_FinishEncode(&rc);
    *out_size = RC_OutSize(&rc);

    free(seq_model);

    return out;
}

char *decode_seq(unsigned char *in,  unsigned int in_size,
		 int ctx_size,
		 /*unsigned char *out,*/ unsigned int out_size) {
    char *out = malloc(out_size);
    if (!out)
	return NULL;

    const int msize = 1<<(2*ctx_size);
    const int mask = msize-1;

    SIMPLE_MODEL(NSYM,_) *seq_model = malloc(msize * sizeof(*seq_model));
    int i;

    // Do histogram to get observed values.
    // Then set m to max number of elements in histogram.
    int m = NSYM;

    for (i = 0; i < msize; i++)
	SIMPLE_MODEL(NSYM,_init)(&seq_model[i], m);

    SIMPLE_MODEL(256,_) run_len4;
    SIMPLE_MODEL(256,_init)(&run_len4, 256);
    SIMPLE_MODEL(256,_) run_lenN;
    SIMPLE_MODEL(256,_init)(&run_lenN, 256);
    SIMPLE_MODEL(256,_) literal;
    SIMPLE_MODEL(256,_init)(&literal, 256);

    RangeCoder rc;
    RC_SetInput(&rc, in, in+in_size);
    RC_StartDecode(&rc);

    int last, last2;
    /* Corresponds to a 12-mer word that doesn't occur in human genome. */
    last  = 0x007616c7 & mask;
    //last2 = (0x2c6b62ff >> (32 - 2*ctx)) & mask; // both strands mode

    // TODO: add an escape mechanism for within ACGT.
    // Like RLE it could be a count of ACGTs, followed by a
    // count of non-ACGTs and a different simpler model.
#if NSYM<=4
    for (int i = 0; i < out_size;) {
	// ACGTs
	int run = 0, r2, j;
	do {
	    r2 = SIMPLE_MODEL(256, _decodeSymbol)(&run_len4, &rc);
	    //fprintf(stderr, "run4 %d\n", r2);
	    run += r2;
	} while (r2 == 255);

	for (j = 0; j < run; j++) {
	    unsigned char b =
		SIMPLE_MODEL(NSYM, _decodeSymbol)(&seq_model[last], &rc);
	    last = ((last<<2) + b) & mask;
	    out[i+j] = "ACGT"[b];
	}
	i += run;
	if (i >= out_size)
	    break;

	// non ACGT
	run = 0;
	do {
	    r2 = SIMPLE_MODEL(256, _decodeSymbol)(&run_lenN, &rc);
	    //fprintf(stderr, "runN %d\n", r2);
	    run += r2;
	} while (r2 == 255);

	for (j = 0; j < run; j++)
	    out[i+j] = SIMPLE_MODEL(256, _decodeSymbol)(&literal, &rc);
	i += run;
    }
#else
    for (int i = 0; i < out_size; i++) {
	unsigned char b =
	    SIMPLE_MODEL(NSYM, _decodeSymbol)(&seq_model[last], &rc);
	last = ((last<<2) + b) & mask;
	out[i] = "ACGTN"[b];
    }
#endif

    RC_FinishDecode(&rc);

    free(seq_model);

    return out;
}


#define BLK_SIZE 300*1000000
//#define BLK_SIZE 10*1000000

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

    // FIXME: use variable sized integers

    // Block based, for arbitrary sizes of input
    if (decomp) {
	for (;;) {
	    // Load next compressed block
	    int i, j, nr, u_len, c_len;
	    char *comp;

	    if (fread(&nr, 1, 4, in_fp) != 4)
		break;
	    fastq *fq = fastq_alloc(nr);

	    // ----------
	    // Name
	    if (fread(&u_len, 1, 4, in_fp) != 4)
		break;
	    if (fread(&c_len, 1, 4, in_fp) != 4)
		break;
	    comp = malloc(c_len);
	    if (fread(comp, 1, c_len, in_fp) != c_len)
		break;

	    out = tok3_decode_names(comp, c_len, &u_len);
	    fq->name_buf = out;
	    fq->name_len = u_len;
	    free(comp);

	    // ----------
	    // Lengths
	    int c, err = 0;
	    if ((c = getc(in_fp)) == EOF) break;
	    if (c > 0) {
		// Fixed length
		char buf[10];
		err |= fread(buf, 1, c, in_fp) != c;
		uint32_t len;
		err |= var_get_u32(buf, buf+c, &len) == 0;
		for (i = 0; i < nr; i++)
		    fq->len[i] = len;
	    } else {
		// Variable length
		uint32_t blen;
		if (fread(&blen, 1, 4, in_fp) != 4)
		    break;
		char *buf = malloc(blen);
		err |= fread(buf, 1, blen, in_fp) != blen;
		int nb = 0;
		for (i = 0; i < nr; i++) {
		    int vl;
		    nb += (vl = var_get_u32(buf+nb, buf+blen, &fq->len[i]));
		    err |= vl == 0;
		    //fread(&fq->len[i], 1, 4, in_fp) != 4;
		}
		if (nb != blen)
		    break;
		free(buf);
	    }
	    if (err)
		break;

	    // ----------
	    // Seq
	    if (fread(&u_len, 1, 4, in_fp) != 4)
		break;
	    if (fread(&c_len, 1, 4, in_fp) != 4)
		break;
	    comp = malloc(c_len);
	    if (fread(comp, 1, c_len, in_fp) != c_len)
		break;
#if 1
	    out = decode_seq(comp, c_len, SEQ_CTX, u_len);
#else
	    out = rans_uncompress_4x16(comp, c_len, &u_len);
#endif
	    free(comp);
	    fq->seq_buf = out;
	    fq->seq_len = u_len;
	    for (i = 0; i < nr; i++)
		fq->seq[i] = i ? fq->seq[i-1] + fq->len[i-1] : 0;

	    // ----------
	    // Qual
	    int mode = getc(in_fp);
	    if (fread(&u_len, 1, 4, in_fp) != 4)
		break;
	    if (fread(&c_len, 1, 4, in_fp) != 4)
		break;
	    comp = malloc(c_len);
	    if (fread(comp, 1, c_len, in_fp) != c_len)
		break;

	    if (mode == 0) {
		// Rans
		out = rans_uncompress_4x16(comp, c_len, &u_len);
		fq->qual_buf = out;
		fq->qual_len = u_len;
	    } else {
		// FQZComp qual
		fqz_slice s;
		s.num_records = fq->num_records;
		s.len = fq->len;
		s.flags = fq->flag;
		//s.seq = (unsigned char **)fq->seq;
		s.seq = (unsigned char **)malloc(fq->num_records * sizeof(char *));
		for (i = j = 0; i < fq->num_records; j += fq->len[i++])
		    s.seq[i] = fq->seq_buf + j;

		int *lengths = malloc(nr * sizeof(lengths));
		out = fqz_decompress((char *)comp, c_len, &out_len,
				     lengths, nr, &s);
		fq->qual_buf = out;
		fq->qual_len = out_len;
		free(s.seq);
		free(lengths);
	    }
	    free(comp);
	    for (i = 0; i < fq->qual_len; i++)
		fq->qual_buf[i] += 33;

	    // ----------
	    // Convert back to fastq
	    char *np = fq->name_buf;
	    char *sp = fq->seq_buf;
	    char *qp = fq->qual_buf;

	    for (i = 0; i < nr; i++) {
		fprintf(out_fp, "@%s\n%.*s\n+\n%.*s\n",
			np, fq->len[i], sp, fq->len[i], qp);
		np += strlen(np)+1;
		sp += fq->len[i];
		qp += fq->len[i];
	    }

	    fastq_free(fq);
	}

    } else {
	// Load the next batch of fastq entries
	for(;;) {
	    fastq *fq = load_seqs(in_fp, blk_size);
	    if (!fq)
		exit(1);

	    // FIXME: report actual block size consumed?
	    fprintf(stderr, "blk_size = %d, nrec = %d\n",
		    blk_size, fq->num_records);
	    if (!fq->num_records) {
		fastq_free(fq);
		break;
	    }

	    //fastq_dump(fq);
	    
	    fwrite(&fq->num_records, 1, 4, out_fp);

	    //----------
	    // Names: tok3
	    int clen;
	    out = tok3_encode_names(fq->name_buf, fq->name_len, 5, 0,
				    &clen, NULL);
	    fwrite(&fq->name_len, 1, 4, out_fp);
	    fwrite(&clen, 1, 4, out_fp);
	    fwrite(out, 1, clen, out_fp);
	    free(out);
	    fprintf(stderr, "Names: %10d to %10d\n", fq->name_len, clen);

	    //----------
	    // Read lengths
	    if (fq->fixed_len) {
		// Fixed length, with next byte holding the size of length
		char buf[5], nb = 1;
		nb += var_put_u32(buf+1, NULL, fq->fixed_len);
		buf[0] = nb-1;
		fwrite(buf, 1, nb, out_fp);
		fprintf(stderr, "Len:   %10d to %10d\n", 4*fq->num_records, nb);
	    } else {
		// Variable length (next byte 0), with 4 byte len followed
		// by var-int lengths.
		int i, nb = 0;
		char *buf = malloc(fq->num_records*5);
		putc(0, out_fp);

		for (i = 0; i < fq->num_records; i++)
		    nb += var_put_u32(buf+nb, NULL, fq->len[i]);
		fwrite(&nb, 1, 4, out_fp);
		fwrite(buf, 1, nb, out_fp);
		free(buf);
		fprintf(stderr, "Len:   %10d to %10d\n",
			4*fq->num_records, nb+5);
	    }

	    //----------
	    // Seq: rans or statistical modelling
#if 1
	    out = encode_seq(fq->seq_buf, fq->seq_len, SEQ_CTX, &clen);
	    fwrite(&fq->seq_len, 1, 4, out_fp);
	    fwrite(&clen, 1, 4, out_fp);
	    fwrite(out, 1, clen, out_fp);
	    free(out);
#else
	    // FIXME: test 197, 65 and 1
	    // FIXME: also try fqz for encoding, qmap and qshift=2.
	    // But need extension to permit reverse complement.
	    // Also not quite the same in terms of model updates?
	    out = rans_compress_4x16(fq->seq_buf, fq->seq_len, &clen, 197);
	    fwrite(&fq->seq_len, 1, 4, out_fp);
	    fwrite(&clen, 1, 4, out_fp);
	    fwrite(out, 1, clen, out_fp);
	    free(out);
#endif
	    fprintf(stderr, "Seq:   %10d to %10d\n", fq->seq_len, clen);

	    //----------
	    // Qual: rans or fqz
	    // Convert fastq struct to fqz_slice for context
	    if (0) {
		// Fast mode
		putc(0, out_fp);
		out = rans_compress_4x16(fq->qual_buf, fq->seq_len,
					 &clen, 193);
		// crash on /tmp/_1m4.fq with mode 197.  Why?
		out_len = clen;
	    } else {
		putc(1, out_fp);
		fqz_slice *s = malloc(fq->num_records * sizeof(*s));
		s->num_records = fq->num_records;
		s->len = fq->len;
		s->flags = fq->flag;
		s->seq = malloc(fq->num_records * sizeof(char *));
		int i, j;
		for (i = j = 0; i < fq->num_records; j += fq->len[i++])
		    s->seq[i] = fq->seq_buf + j;

		// FIXME: expose fqz_pick_parameters function so we
		// can initialise it here and then also turn off DO_LEN.

		// Concatenate qualities together into a single block.
		// FIXME: move qual down by 33, '!'.
		out = fqz_compress(vers, s, fq->qual_buf, fq->qual_len,
				   &out_len, strat, gp);
		free(s->seq);
		free(s);
	    }
	    fwrite(&fq->qual_len, 1, 4, out_fp);
	    fwrite(&out_len, 1, 4, out_fp);
	    fwrite(out, 1, out_len, out_fp);
	    free(out);
	    fprintf(stderr, "Qual:  %10d to %10d\n", fq->qual_len, (int)out_len);

	    fastq_free(fq);
	}
    }

    fclose(in_fp);
    fclose(out_fp);

    return 0;
}
