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

/*
File format:

[Header]  TODO (magic number)
[Block]*  Zero or more blocks of records
[Index]   TODO: self index + metrics

Where Block is:

4    Num records
4    TODO: Total compressed block size
[Block CRC]  TODO

1    Name strategy
4    NU: uncompressed name length
4    NC: compressed name length
NC   Name compressed data

1    Read length strategy
?    Read length data

1    Sequence strategy (bits 0..2), both_strands (bit 3), level (4..7)
4    SU: Uncompressed sequence size
4    SC: Compressed sequence size
SC   Compressed sequence data

1    Quality strategy
4    QU: Uncompressed quality size
4    QC: Compressed quality size
QC   Compressed quality data


 */

// TODO
// - Split aux tags into own data series using CRAM's TL + per tag.

// - Seq encoding using STR + copy-number?
//   Couldn't get this to work well though, even though it sounds like an
//   ideal thing for ONT/PB-CLR.  Maybe let it kick in after so many bases
//   in a homopolymer?

// - Removal of lengths from fqzqual stream

// - Entropy encoding of read length stream

// - Exploration of multiple rans options (o4, o5, o133, o197?)

// - Reuse of memory buffers for speed

// - Increase sizes of fqzcomp contexts?  Not so useful for CRAM, but maybe
//   it's still beneficial to go beyond 16-bit.

// - Improve fqzcomp tables to permit any mapping rather than monotonic.

// - Distinguish explicit method opts (-s1 -S13B -q1 -Q2 etc) from auto
//   picked options (-3, -5) which use auto-selected metrics

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <fcntl.h>
#include <ctype.h>
#include <limits.h>
#include <sys/time.h>

#include "htscodecs/varint.h"
#include "htscodecs/fqzcomp_qual.h"
#include "htscodecs/tokenise_name3.h"
#include "htscodecs/rANS_static4x16.h"
#include "htscodecs/varint.h"
#include "thread_pool.h"
#include "lzp16e.h"

#define BLK_SIZE 512*1000000

#ifndef MIN
#  define MIN(a,b) ((a)<(b)?(a):(b))
#  define MAX(a,b) ((a)>(b)?(a):(b))
#endif

// Review metrics stats every X blocks for Y trials
#define METRICS_REVIEW 50
#define METRICS_TRIAL 3

typedef enum {
    SEC_NAME,
    SEC_LEN,
    SEC_SEQ,
    SEC_QUAL,
    SEC_LAST
} sections;

typedef enum {
    // general purpose
    RANS0=1, RANS1, RANS64, RANS65, RANS128, RANS129, RANS192, RANS193,
    //RANSXN1  TODO: stripe N-way where N is fixed read len?

    // LZP; differing min lengths?  Make len part of format?
    LZP3, /* TODO: could use on seq too maybe */
    // TODO LZP2, LZP4, LZP16? Needs storing in byte stream too

    // Name specific
    TOK3_3, TOK3_5, TOK3_7, TOK3_9,
    TOK3_3_LZP, TOK3_5_LZP, TOK3_7_LZP, TOK3_9_LZP,

    // Seq
    SEQ10, SEQ12, SEQ12B, SEQ13B, SEQ14B, SEQ_CUSTOM/*TODO*/,

    // Qual
    FQZ_AUTO, FQZ0, FQZ1, FQZ2, FQZ3, QUAL_CUSTOM/*TODO*/,

    M_LAST,
} methods;

typedef struct {
    uint64_t usize[M_LAST], csize[M_LAST];   // current accumulated sizes
    int review, trial;
} metrics;

static uint32_t method_avail[SEC_LAST];
static methods method_used[SEC_LAST];
static metrics stats[SEC_LAST];

typedef struct {
    int num_records;                 // number of fastq entries
    char *name_buf;                  // concatenated names, \0 separator
    char *seq_buf;                   // concatenated seq,   no separator
    char *qual_buf;                  // concatenated qual,  no separator
    int *name;                       // index into name_buf
    int *seq;                        // index into seq_buf
    int *qual;                       // index into qual_buf
    unsigned int *len;               // sequence length
    unsigned int *flag;              // READ1/READ2 parsed from name
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

#define goto if (fprintf(stderr, "ERR %s:%d\n", __FILE__, __LINE__)) goto
fastq *load_seqs(char *in, int blk_size, int *last_offset) {
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
    int last_start = 0;
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
	c = in[i++];
	if (c != '@')
	    goto err;

	int name_i_ = name_i; // tmp copy so we can unwind a partial decode
	while (i < blk_size && (c = in[i++]) && c != '\n') {
	    if (name_i_ >= name_sz) {
		name_sz = name_sz * 1.5 + 1000;
		name_buf = fq->name_buf = realloc(fq->name_buf, name_sz);
	    }
	    name_buf[name_i_++] = c;
	}
	if (i == blk_size)
	    break;

	name_buf[name_i_++] = 0;

	int flag = 0;
	if (name_i_ > 3 &&
	    name_buf[name_i_-2] == '2' &&
	    name_buf[name_i_-3] == '/')
	    flag = FQZ_FREAD2;
	if (last_name >= 0 &&
	    strcmp(fq->name_buf + fq->name[nr], fq->name_buf + last_name) == 0)
	    flag = FQZ_FREAD2;
	fq->flag[nr] = flag;
	last_name = fq->name[nr];

	// seq
	fq->seq[nr] = seq_i;
	int len = seq_i;
	int seq_i_ = seq_i;
	while (i < blk_size && (c = in[i++]) && c != '\n') {
	    if (seq_i_ >= seq_sz) {
		// very unlikely given blk_size/2 starting point,
		// but not impossible.
		seq_sz = seq_sz * 1.5 + 1000;
		seq_buf = fq->seq_buf = realloc(fq->seq_buf, seq_sz);
	    }
	    seq_buf[seq_i_++] = c;
	}
	if (i == blk_size)
	    break;
	fq->len[nr] = seq_i_ - len;
	//seq_buf[seq_i_++] = '\n'; // use instead of length terminator?

	if (fq->fixed_len == -1)
	    fq->fixed_len = fq->len[nr];
	else if (fq->fixed_len > 0)
	    if (fq->fixed_len != fq->len[nr])
		fq->fixed_len = 0;

	// +(name)
	if (i < blk_size && (c = in[i++]) != '+')
	    goto err;
	while (i < blk_size && (c = in[i++]) && c != '\n')
	    ;
	if (i == blk_size)
	    break;

	// qual
	fq->qual[nr] = qual_i;
	len = qual_i;
	int qual_i_ = qual_i;
	while (i < blk_size && (c = in[i++]) && c != '\n') {
	    if (qual_i_ >= qual_sz) {
		// very unlikely given blk_size/2 starting point
		qual_sz = qual_sz * 1.5 + 1000;
		qual_buf = fq->qual_buf = realloc(fq->qual_buf, qual_sz);
	    }
	    qual_buf[qual_i_++] = c-33;
	}

	if (fq->len[nr] != qual_i_ - len) {
	    if (i == blk_size)
		break;

	    goto err;
	} else if (i == blk_size && c != '\n')
	    break;

	name_i = name_i_;
	seq_i  = seq_i_;
	qual_i = qual_i_;

	last_start = i;
	nr++;
    }

    // FIXME: need to deal with case where block size is smaller than
    // a single record!  For now that puts a limit on smallest block size.
    *last_offset = last_start; // so we can continue for next block
    
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

#if 0
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
#endif

#define NSYM 256
#define STEP 8
#include "htscodecs/c_simple_model.h"
#undef NSYM
#undef STEP

#define NSYM 4
#include "htscodecs/c_small_model.h"

#undef NSYM
#define NSYM 2
#include "htscodecs/c_small_model.h"

// An order-N arithmetic encoder, dedicated to sequence contexts.
char *encode_seq(unsigned char *in,  unsigned int in_size,
		 unsigned int *len, int nrecords, int both_strands,
		 int ctx_size, unsigned int *out_size) {
    char *out = malloc(in_size + 100);
    if (!out)
	return NULL;

    const int msize = 1<<(2*ctx_size);
    const int mask = msize-1;
    int i;

    SMALL_MODEL(4,_) *seq_model = malloc(msize * sizeof(*seq_model));

    for (i = 0; i < msize; i++)
	SMALL_MODEL(4,_init)(&seq_model[i]);

    SMALL_MODEL(2,_) state_model[3];
    SIMPLE_MODEL(256,_) run_len[3];
    for (i = 0; i < 3; i++) {
	SMALL_MODEL(2,_init)(&state_model[i]);
	SIMPLE_MODEL(256,_init)(&run_len[i], 256);
    }

    SIMPLE_MODEL(256,_) literal;
    SIMPLE_MODEL(256,_init)(&literal, 256);

    RangeCoder rc;
    RC_SetOutput(&rc, (char *)out);
    RC_StartEncode(&rc);

    /* Corresponds to a 12-mer word that doesn't occur in human genome. */
    int last  = 0x007616c7 & mask;
    int last2 = (0x2c6b62ff >> (32 - 2*ctx_size)) & mask; // both strands mode

    int L[256];
    for (int i = 0; i < 256; i++)
	L[i] = 4; // N
    L['A'] = 0;
    L['C'] = 1;
    L['G'] = 2;
    L['T'] = 3;

    L['a'] = 0x80;
    L['c'] = 0x81;
    L['g'] = 0x82;
    L['t'] = 0x83;

    // Transition table to stored code:
    //    uc lc N
    // uc -  0  1
    // lc 0  -  1
    // N  0  1  -
    enum { uc_ACGT = 0, lc_ACGT = 1, other = 2 } state = uc_ACGT;

    int nseq = 0;
    int seq_len = len[nseq++];
    for (int i = 0; i < in_size; ) {//i++) {
	// Count size of consecutive symbols matching the same state
	int j, run, r2;
	switch (state) {
	case uc_ACGT: // uppercase ACGT symbols
	    for (j = i; j < in_size; j++) {
		if (L[in[j]] >= 4)
		    break;
	    }
	    break;

	case lc_ACGT: // lowercase acgt symbols
	    for (j = i; j < in_size; j++) {
		if (L[in[j]] < 0x80)
		    break;
	    }
	    break;

	case other: // ambiguity codes
	    for (j = i; j < in_size; j++) {
		if (L[in[j]] != 4)
		    break;
	    }
	    break;
	}

	// Encode the run length
	r2 = run = j-i;
	for (;;) {
	    //fprintf(stderr, "Encode %d of %d for state %d\n",
	    //    MIN(255, r2), run, state);
	    SIMPLE_MODEL(256, _encodeSymbol)(&run_len[state], &rc,
					     MIN(255, r2));
	    if (r2 >= 255)
		r2 -= 255;
	    else
		break;
	}

	// Encode the symbols
	switch (state) {
	case uc_ACGT:
	case lc_ACGT:
	    for (j = 0; j < run; j++) {
		unsigned char b = L[in[i+j]] & 3;
		SMALL_MODEL(4, _encodeSymbol)(&seq_model[last], &rc, b);

		last = ((last<<2) + b) & mask;
		int pf = ((last<<6)&mask)
			 +(i+j+3<in_size
			   ?L[in[i+j+1]]*16+L[in[i+j+2]]*4+L[in[i+j+3]]
			   :0);
		_mm_prefetch((const char *)&seq_model[pf], _MM_HINT_T0);

		// 0.7% and 3.2% smaller for _.FQ and _.fq respectively
		// (at ctx_size 12), but 45% more CPU for seq encoding.
		if (both_strands) {
		    int b2 = last2 & 3;
		    last2 = last2/4 + ((3-b) << (2*ctx_size-2));
		    SMALL_MODEL(4, _updateSymbol)(&seq_model[last2], b2);

		    // ~25% speed gain by prefetching bottom strand too
		    uint32_t i3 = i+j+3 < in_size
			? L[in[i+j+1]] + L[in[i+j+2]]*4 + L[in[i+j+3]]*16
			: 0;
		    i3 = (0x3f - i3) << (2*ctx_size-6);
		    pf = i+j+3 < in_size ? (last2>>6) +i3 : 0;
		    _mm_prefetch((const char *)&seq_model[pf], _MM_HINT_T0);
		}

		// In theory we should reset context for each new sequence
		// as there is no obvious correlation between one sequence
		// and the next.  In practice the difference is 1-2%.
		// It only costs 1-2% CPU too, so worth doing.
		//
		//          -s5/100MB             -s7 -b / 512MB
		// Without: 76655812              60731188
		// With:    75789505 -1.1%        59638082 -1.8%
		// Slowdown 0.2%                  1.9%
		if (--seq_len == 0 && i+j+1 < in_size) {
		    if (nseq >= nrecords) {
			free(out);
			free(seq_model);
			return NULL;
		    }
		    seq_len = len[nseq++];
		    last = 0x007616c7 & mask;
		    last2 = (0x2c6b62ff >> (32 - 2*ctx_size)) & mask;
		}
	    }
	    break;

	case other:
	    for (j = 0; j < run; j++) {
		SIMPLE_MODEL(256, _encodeSymbol)(&literal, &rc, in[i+j]);
		if (--seq_len == 0 && i+j+1 < in_size) {
		    if (nseq >= nrecords) {
			free(out);
			free(seq_model);
			return NULL;
		    }
		    seq_len = len[nseq++];
		    last = 0x007616c7 & mask;
		    last2 = (0x2c6b62ff >> (32 - 2*ctx_size)) & mask;
		}
	    }
	}

	i += run;
	if (i >= in_size)
	    break;

	// Encode switch to next state
	switch(L[in[i]]) {
	case 0: case 1: case 2: case 3:
	    //fprintf(stderr, "state %d to 0, => 0\n", state);
	    SMALL_MODEL(2, _encodeSymbol)(&state_model[state], &rc, 0);
	    state = uc_ACGT;
	    break;
	case 0x80: case 0x81: case 0x82: case 0x83:
	    //fprintf(stderr, "state %d to 1, => %d\n", state, state==other);
	    SMALL_MODEL(2, _encodeSymbol)(&state_model[state], &rc,
					  state == other);
	    state = lc_ACGT;
	    break;
	default:
	    //fprintf(stderr, "state %d to 2, => 1\n", state);
	    SMALL_MODEL(2, _encodeSymbol)(&state_model[state], &rc, 1);
	    state = other;
	    break;
	}
    }

    RC_FinishEncode(&rc);
    *out_size = RC_OutSize(&rc);

    free(seq_model);

    return out;
}

char *decode_seq(unsigned char *in,  unsigned int in_size,
		 unsigned int *len, int nrecords, int both_strands,
		 int ctx_size, unsigned int out_size) {
    char *out = malloc(out_size);
    if (!out)
	return NULL;

    const int msize = 1<<(2*ctx_size);
    const int mask = msize-1;
    int i;

    SMALL_MODEL(4,_) *seq_model = malloc(msize * sizeof(*seq_model));

    // Do histogram to get observed values.
    // Then set m to max number of elements in histogram.
    for (i = 0; i < msize; i++)
	SMALL_MODEL(4,_init)(&seq_model[i]);

    SMALL_MODEL(2,_) state_model[3];
    SIMPLE_MODEL(256,_) run_len[3];
    for (i = 0; i < 3; i++) {
	SMALL_MODEL(2,_init)(&state_model[i]);
	SIMPLE_MODEL(256,_init)(&run_len[i], 256);
    }

    SIMPLE_MODEL(256,_) literal;
    SIMPLE_MODEL(256,_init)(&literal, 256);

    RangeCoder rc;
    RC_SetInput(&rc, (char *)in, (char *)in+in_size);
    RC_StartDecode(&rc);

    /* Corresponds to a 12-mer word that doesn't occur in human genome. */
    int last  = 0x007616c7 & mask;
    int last2 = (0x2c6b62ff >> (32 - 2*ctx_size)) & mask; // both strands mode

    // Transition table to stored code:
    //    uc lc N
    // uc -  0  1
    // lc 0  -  1
    // N  0  1  -
    enum { uc_ACGT = 0, lc_ACGT = 1, other = 2 } state = uc_ACGT;

    int nseq = 0;
    int seq_len = len[nseq++];
    for (int i = 0; i < out_size;) {
	int j, run = 0, r2;
	// Fetch run length
	do {
	    r2 = SIMPLE_MODEL(256, _decodeSymbol)(&run_len[state], &rc);
	    run += r2;
	    //fprintf(stderr, "Decode %d of %d for state %d\n", r2, run, state);
	} while (r2 == 255);

	if (i + run > out_size)
	    // or error as it's malformed data
	    run = out_size - i;

	// Decode
	switch (state) {
	case uc_ACGT:
	case lc_ACGT: {
	    char *bases = state==lc_ACGT ? "acgt" : "ACGT";
	    for (j = 0; j < run; j++) {
		unsigned char b =
		    SMALL_MODEL(4, _decodeSymbol)(&seq_model[last], &rc);
		last = ((last<<2) + b) & mask;
		_mm_prefetch((const char *)&seq_model[(last<<4)&mask],
			     _MM_HINT_T0);
		out[i+j] = bases[b];

		if (both_strands) {
		    int b2 = last2 & 3;
		    last2 = last2/4 + ((3-b) << (2*ctx_size-2));
		    SMALL_MODEL(4, _updateSymbol)(&seq_model[last2], b2);
	        }

		if (--seq_len == 0 && i+j+1 < out_size) {
		    if (nseq >=  nrecords) {
			free(out);
			free(seq_model);
			return NULL;
		    }
		    seq_len = len[nseq++];
		    last = 0x007616c7 & mask;
		    last2 = (0x2c6b62ff >> (32 - 2*ctx_size)) & mask;
		}
	    }
	    break;
	}

	case other: // ambiguity codes
	    for (j = 0; j < run; j++) {
		out[i+j] = SIMPLE_MODEL(256, _decodeSymbol)(&literal, &rc);
		if (--seq_len == 0 && i+j+1 < out_size) {
		    if (nseq >=  nrecords) {
			free(out);
			free(seq_model);
			return NULL;
		    }
		    seq_len = len[nseq++];
		    last = 0x007616c7 & mask;
		    last2 = (0x2c6b62ff >> (32 - 2*ctx_size)) & mask;
		}
	    }
	    break;
	}

	i += run;
	if (i >= out_size)
	    break;

	// Next state
	int nstate = SMALL_MODEL(2, _decodeSymbol)(&state_model[state], &rc);
	switch (state) {
	case uc_ACGT:
	    state = nstate ? other : lc_ACGT;
	    break;
	case lc_ACGT:
	    state = nstate ? other : uc_ACGT;
	    break;
	case other:
	    state = nstate ? lc_ACGT : uc_ACGT;
	    break;
	}
    }

    RC_FinishDecode(&rc);

    free(seq_model);

    return out;
}

static char *encode_names(unsigned char *name_buf,  unsigned int name_len,
			  int strat, int level, unsigned int *out_size) {
    // TODO: work out a better maximum bound
    char *nout = malloc(name_len*2+1000), *cp = nout;
    
    *(uint32_t *)cp = name_len; cp += 4;
    *cp++ = strat; // name method;

    if (strat == 0) {
	unsigned char *lzp_out = (unsigned char *)cp;
	unsigned int clen = lzp(name_buf, name_len, lzp_out);
	unsigned char *out = rans_compress_4x16(lzp_out, clen, &clen, 5);
	*(uint32_t *)cp = clen; cp += 4;
	memcpy(cp, out, clen);  cp += clen;
	free(out);

    } else if (strat == 1) {
	int clen;
	unsigned char *out = tok3_encode_names((char *)name_buf, name_len,
					       level, 0, &clen, NULL);
	*(uint32_t *)cp = clen; cp += 4;
	memcpy(cp, out, clen);  cp += clen;
	free(out);

    } else {
	char *n1 = malloc(name_len);
	char *n2 = malloc(name_len);
	unsigned char *flag = malloc(name_len/2);  //Worst case\n 
	char *cp1 = n1, *cp2 = n2;
	int i = 0, nr = 0;
	// Flag bit 0: has "/NUM"
	// Flag bit 1: /1 vs /2
	// Flag bit 2: has a comment
	// Flag bit 3: space vs tab before comment
	while (i < name_len) {
	    int j, f = 0;
	    int w1end = 0;
	    int w2start = 0;
	    int w2end = 0;
	    for (j = i; j < name_len; j++) {
		if (name_buf[j] == '\0') {
		    w2end = j;
		    break;
		}
		if (!w2start && (name_buf[j] == ' ' ||
				 name_buf[j] == '\t')) {
		    w1end = j;
		    w2start = j+1;
		    f |= 4; // FLAG: has comment
		}
	    }

	    if (!w1end)
		w1end = j;
	    if (!w2end)
		w2end = j;

	    if (w2start)
		// FLAG: space vs tab
		f |= name_buf[w2start-1] == ' ' ? 0 : 8;

	    if (w1end>1 && name_buf[w1end-2] == '/') {
		// FLAG /1 or /2
		if (name_buf[w1end-1] == '1')
		    f |= 1, w1end -= 2;
		else if (name_buf[w1end-1] == '2')
		    f |= 3, w1end -= 2;
	    }

	    flag[nr++] = f;
	    memcpy(cp1, &name_buf[i], w1end-i);
	    cp1[w1end-i]=0;
	    cp1 += w1end-i+1;

	    if (w2start) {
		memcpy(cp2, &name_buf[w2start], w2end-w2start);
		cp2[w2end-w2start] = 0;
		cp2 += w2end-w2start+1;
	    }

	    i = j+1;
	}

	int clen1;
	unsigned int clenf, clen2 = 0;
	unsigned char *out = tok3_encode_names(n1, cp1-n1, level, 0, &clen1,
					       NULL);
	unsigned char *outf = rans_compress_4x16(flag, nr, &clenf, 129);
	unsigned char *out2 = NULL;
	if (cp2 != n2) {
	    unsigned char *lzp_out = malloc((cp2-n2)*2);
	    clen2 = lzp((unsigned char *)n2, cp2-n2, lzp_out);
	    out2 = rans_compress_4x16(lzp_out, clen2, &clen2, 5);
	    free(lzp_out);
	}

	unsigned int clen = clen1 + clenf + clen2 + 8;

	*(uint32_t *)cp = clen;  cp += 4;
	*(uint32_t *)cp = clen1; cp += 4;
	*(uint32_t *)cp = clenf; cp += 4;

	memcpy(cp, out, clen1);  cp += clen1;
	free(out);

	memcpy(cp, outf, clenf); cp += clenf;
	free(outf);

	if (out2) {
	    memcpy(cp, out2, clen2); cp += clen2;
	    free(out2);
	}
	free(n1);
	free(n2);
	free(flag);
    }

    *out_size = cp - nout;
    return nout;
}

static char *decode_names(unsigned char *comp,  unsigned int c_len,
			  unsigned int u_len, int strat) {
    unsigned char *out;

    if (strat == 0) {
	unsigned int ru_len;
	unsigned char *rout = rans_uncompress_4x16(comp, c_len, &ru_len);
	out = malloc(u_len);
	u_len = unlzp(rout, ru_len, out);
	free(rout);
    } else if (strat == 1) {
	out = tok3_decode_names(comp, c_len, &u_len);
    } else {
	uint32_t clen1 = *(uint32_t *)comp;
	uint32_t clenf = *(uint32_t *)(comp+4);
	uint32_t clen2 = c_len - clen1 - clenf - 8;

	// Uncompress 3 separate components
	unsigned int u_len1, u_lenf, u_len2;
	unsigned char *out1 = tok3_decode_names(comp+8, clen1, &u_len1);
	unsigned char *outf = rans_uncompress_4x16(comp+8+clen1, clenf,
						   &u_lenf);
	unsigned char *out2 = NULL;
	if (clen2) {
	    unsigned int rulen;
	    unsigned char *rout = rans_uncompress_4x16(comp+8+clen1+clenf,
						       clen2, &rulen);
	    out2 = malloc(u_len);
	    u_len2 = unlzp(rout, rulen, out2);
	    free(rout);
	}

	// Stitch together ID + flag + comment
	unsigned char *cp1 = out1, *cp1_end = out1+u_len1;
	unsigned char *cpf = outf, *cpf_end = outf+u_lenf;
	unsigned char *cp2 = out2, *cp2_end = out2 + u_len2;
	out = malloc(u_len);
	unsigned char *cp  = out,  *cp_end = out + u_len;
	unsigned char *last_cp = NULL;
	while (cp < cp_end) {
	    while (cp1 < cp1_end && cp < cp_end && *cp1)
		*cp++ = *cp1++;
	    cp1++;

	    int flag = 0;
	    if (cpf < cpf_end)
		flag = *cpf++;
	    if ((flag & 1) && cp+1 < cp_end) {
		*cp++ = '/';
		*cp++ = (flag & 2) ? '2' : '1';
	    }
		
	    if (flag & 4)
		*cp++ = (flag & 8) ? '\t' : ' ';

	    if (cp2) {
		while (cp2 < cp2_end && cp < cp_end && *cp2)
		    *cp++ = *cp2++;
		cp2++;
	    }

	    if (cp == last_cp)
		// ran out of data early; avoids looping forever
		break;

	    if (cp < cp_end) {
		*cp++ = 0;
	    } else {
		free(out1);
		free(outf);
		free(out2);
		goto err;
	    }
	    last_cp = cp;
	}

	free(out1);
	free(outf);
	free(out2);
    }

    return (char *)out;

 err:
    free(out);
    return NULL;
}

typedef struct {
    int nstrat, sstrat, qstrat;
    int nlevel, slevel, qlevel;
    int verbose;
    int both_strands;
    uint32_t blk_size;
    int nthread;
} opts;

typedef struct {
    int64_t nblock;
    int64_t nusize, ncsize, ntime;
    int64_t susize, scsize, stime;
    int64_t qusize, qcsize, qtime;
    int64_t lusize, lcsize, ltime;
} timings;

static inline uint64_t tvdiff(struct timeval *tv1, struct timeval *tv2) {
    return (tv2->tv_sec - tv1->tv_sec) * 1000000
	+ tv2->tv_usec - tv1->tv_usec;
}

void update_stats(timings *t,
		  int column, // 0=name 1=seq 2=qual 3=length
		  int64_t usize, int csize, int time) {
    switch (column) {
    case 0:
	t->nusize += usize;
	t->ncsize += csize;
	t->ntime  += time;
	break;
    case 1:
	t->susize += usize;
	t->scsize += csize;
	t->stime  += time;
	break;
    case 2:
	t->qusize += usize;
	t->qcsize += csize;
	t->qtime  += time;
	break;
    case 3:
	t->lusize += usize;
	t->lcsize += csize;
	t->ltime  += time;
	break;
    }
}

void append_timings(timings *t1, timings *t2, int verbose) {
    t1->nblock++;

    t1->nusize += t2->nusize;
    t1->ncsize += t2->ncsize;
    t1->ntime  += t2->ntime;

    t1->susize += t2->susize;
    t1->scsize += t2->scsize;
    t1->stime  += t2->stime;

    t1->qusize += t2->qusize;
    t1->qcsize += t2->qcsize;
    t1->qtime  += t2->qtime;

    t1->lusize += t2->lusize;
    t1->lcsize += t2->lcsize;
    t1->ltime  += t2->ltime;

    if (verbose) {
	fprintf(stderr, "Names   %10ld to %10ld in %.2f sec\n",
		t2->nusize, t2->ncsize, t2->ntime/1e6);
	fprintf(stderr, "Lengths %10ld to %10ld in %.2f sec\n",
		t2->lusize, t2->lcsize, t2->ltime/1e6);
	fprintf(stderr, "Seqs    %10ld to %10ld in %.2f sec\n",
		t2->susize, t2->scsize, t2->stime/1e6);
	fprintf(stderr, "Quals   %10ld to %10ld in %.2f sec\n\n",
		t2->qusize, t2->qcsize, t2->qtime/1e6);
    }
}

#define APPEND_OUT(dat, len) do {			\
    if (*out_size < comp_sz + (len)) {			\
        *out_size = (comp_sz + (len))*1.5 + 1000;	\
	comp = realloc(comp, *out_size);		\
    }							\
    memcpy(comp+comp_sz, (dat), (len));			\
    comp_sz += (len);					\
} while(0);


// Updates the metrics counters and returns the methods to use.
// Method returned is a bitfield of (1<<method_num).
int metrics_method(int sec) {
    // FIXME: add locking
    if (stats[sec].review == 0) {
	stats[sec].review = METRICS_REVIEW;
	stats[sec].trial  = METRICS_TRIAL;
	memset(stats[sec].usize, 0, M_LAST * sizeof(*stats[sec].usize));
	memset(stats[sec].csize, 0, M_LAST * sizeof(*stats[sec].csize));
    }

    int method;
    if (stats[sec].trial>0) {
	method = method_avail[sec];
	stats[sec].trial--;
    } else if (stats[sec].trial == 0) {
	stats[sec].trial--;
	int m, best_m = 0;
	uint32_t best_sz = UINT_MAX;
	for (m = 0; m < M_LAST; m++) {
	    // TODO: parameterise by speed too, plus small block offset?
	    if (best_sz > stats[sec].csize[m] && stats[sec].csize[m]) {
		best_sz = stats[sec].csize[m];
		best_m = m;
	    }
	}
	fprintf(stderr, "Choose best method %d for sec %d\n", best_m, sec);
	method_used[sec] = best_m;
	method = 1<<method_used[sec];

	// TODO: methods that consistently get rejected can be removed from
	// the method_avail.  This is a second level based on
	// total accumulated size stats.
    } else {
	method = 1<<method_used[sec];
    }
    stats[sec].review--;

    return method;
}

// Update the metrics for a given section and method.
// Method parameter isn't a bitfield, but the method number itself.
void metrics_update(int sec, int method, int64_t usize, int64_t csize) {
    stats[sec].usize[method] += usize;
    stats[sec].csize[method] += csize;
    //fprintf(stderr, "Section %d  method %d  size %ld\n", sec, method, csize);
}

// TODO: return buffer + meta, so we don't need memcpy and memmove calls.
char *compress_with_methods(fqz_gparams *gp, fastq *fq, uint32_t methods,
			    int sec, char *in, unsigned int in_size,
			    unsigned int *out_size, int *strat) {
    uint8_t *best_comp = NULL;
    uint32_t best_sz = UINT_MAX;
    int      best_strat = 0;
    char *out;
    int m;
    size_t out_len;

    for (m = 0; m < M_LAST; m++) {
	if (!(methods & (1<<m)))
	    continue;

	switch (m) {
	case RANS0:
	case RANS1:
	case RANS64:
	case RANS65:
	case RANS128:
	case RANS129:
	case RANS192:
	case RANS193: {
	    int order[] = {0,1,64,65,128,129,192,193};
	    *strat = 0;
	    out = (char *)rans_compress_4x16((uint8_t *)in, in_size,
					     out_size, order[m-RANS0]);
	    out_len = *out_size;
	    break;
	}

	case LZP3:
	    out = encode_names(in, in_size, 0 /* LZP + rANS o5 */,
			       (m-TOK3_3)*2+3, out_size);
	    out_len = *out_size;
	    break;

	case TOK3_3:
	case TOK3_5:
	case TOK3_7:
	case TOK3_9:
	    out = encode_names(in, in_size, 1 /* TOK3 */,
			       (m-TOK3_3)*2+3, out_size);
	    out_len = *out_size;
	    break;

	case TOK3_3_LZP:
	case TOK3_5_LZP:
	case TOK3_7_LZP:
	case TOK3_9_LZP:
	    out = encode_names(in, in_size, 2 /* TOK3+LZP */,
			       (m-TOK3_3_LZP)*2+3, out_size);
	    out_len = *out_size;
	    break;

	case SEQ10:
	case SEQ12:
	case SEQ12B:
	case SEQ13B:
	case SEQ14B: {
	    int slevel[]   = {10,12,12,13,14};
	    int both_str[] = {0, 0, 1, 1, 1};

	    int s = m-SEQ10;

	    *strat = (slevel[s]<<4) | (both_str[s]<<3) | 1;
	    out = encode_seq(in, in_size, fq->len, fq->num_records,
			     both_str[s], slevel[s], out_size);
	    out_len = *out_size;
	    break;
	}

	case FQZ0:
	case FQZ1:
	case FQZ2:
	case FQZ3: {
	    *strat = 1;
	    fqz_slice *s = malloc(fq->num_records * sizeof(*s));
	    s->num_records = fq->num_records;
	    s->len = fq->len;
	    s->flags = fq->flag;
	    s->seq = malloc(fq->num_records * sizeof(char *));
	    int i, j;
	    for (i = j = 0; i < fq->num_records; j += fq->len[i++])
		s->seq[i] = (unsigned char *)fq->seq_buf + j;

	    // FIXME: expose fqz_pick_parameters function so we
	    // can initialise it here and then also turn off DO_LEN.
	    out = fqz_compress(4, s, in, in_size, 
			       &out_len, m-FQZ0, gp);
	    free(s->seq);
	    free(s);
	    break;
	}

	default:
	    fprintf(stderr, "Unsupported method %d (set 0x%x)\n", m, methods);
	    abort();
	}

	if (best_sz > out_len) {
	    best_sz = out_len;
	    best_comp = out;
	    best_strat = *strat;
	} else {
	    free(out);
	}
	metrics_update(sec, m, in_size, out_len);
    }
    out = best_comp;

    *out_size = best_sz;
    *strat    = best_strat;
    return out;
}

// Encodes a single block of data
char *encode_block(fqz_gparams *gp, opts *arg, fastq *fq, timings *t,
		   unsigned int *out_size) {
    struct timeval tv1, tv2;

    // Starting guess
    *out_size = 1000;//arg->blk_size/4 + 10000;
    char *comp = malloc(*out_size), *out;
    unsigned int clen, comp_sz = 0;
    int strat = 0, method;

    APPEND_OUT(&fq->num_records, 4);

    //----------
    // Names: tok3
    // Strat 0 = LZP + rANS
    // Strat 1 = Tok3
    // Strat 2 = Name(tok3)+Flag(RC)+Comment(LZP+rANS)
    gettimeofday(&tv1, NULL);
#if 1
    method = metrics_method(SEC_NAME);
    out = compress_with_methods(gp, fq, method, SEC_NAME,
				fq->name_buf, fq->name_len,
				&clen, &strat);
#else
    out = encode_names((uint8_t *)fq->name_buf, fq->name_len, arg->nstrat,
		       arg->nlevel, &clen);
#endif
    APPEND_OUT(out, clen);
    free(out);

    gettimeofday(&tv2, NULL);
    update_stats(t, 0, fq->name_len, clen, tvdiff(&tv1, &tv2));

    //----------
    // Read lengths
    if (fq->fixed_len) {
	// Fixed length, with next byte holding the size of length
	unsigned char buf[5], nb = 1;
	nb += var_put_u32(buf+1, NULL, fq->fixed_len);
	buf[0] = nb-1;
	update_stats(t, 3, 4*fq->num_records, nb, 0);
	APPEND_OUT(buf, nb);
    } else {
	// Variable length (next byte 0), with 4 byte len followed
	// by var-int lengths.
	int i, nb = 0;
	unsigned char *buf = malloc(fq->num_records*5+5);
	buf[nb++] = 0;
	nb += 4; // comp.size placeholder

	for (i = 0; i < fq->num_records; i++)
	    nb += var_put_u32(buf+nb, NULL, fq->len[i]);
	*(uint32_t *)(buf+1) = nb-5;
	APPEND_OUT(buf, nb);
	free(buf);

	update_stats(t, 3, 4*fq->num_records, nb, 0);
    }

    //----------
    // Seq: rans or statistical modelling
    gettimeofday(&tv1, NULL);
    uint8_t  meta[9];

    method = metrics_method(SEC_SEQ);
    out = compress_with_methods(gp, fq, method, SEC_SEQ,
				fq->seq_buf, fq->seq_len,
				&clen, &strat);
    meta[0] = strat;

    *(uint32_t *)(&meta[1]) = fq->seq_len;
    *(uint32_t *)(&meta[5]) = clen;
    APPEND_OUT(meta, 9);
    APPEND_OUT(out, clen);
    free(out);

    gettimeofday(&tv2, NULL);
    update_stats(t, 1, fq->seq_len, clen+9, tvdiff(&tv1, &tv2));

    //----------
    // Qual: rans or fqz
    // Convert fastq struct to fqz_slice for context
    gettimeofday(&tv1, NULL);
    size_t out_len = 0;

    method = metrics_method(SEC_QUAL);
    out = compress_with_methods(gp, fq, method, SEC_QUAL,
				fq->qual_buf, fq->qual_len,
				&clen, &strat);
    meta[0] = strat;
    //fprintf(stderr, "Qual %d -> %d via %d\n", fq->qual_len, clen, strat);

    *(uint32_t *)(&meta[1]) = fq->qual_len;
    *(uint32_t *)(&meta[5]) = clen;
    APPEND_OUT(meta, 9);
    APPEND_OUT(out, clen);
    free(out);

    gettimeofday(&tv2, NULL);
    update_stats(t, 2, fq->qual_len, clen+9, tvdiff(&tv1, &tv2));

    *out_size = comp_sz;

    return comp;
}

#define GET(ptr, len)				\
    do {					\
        if (in_off + (len) > in_size)		\
	    goto err;				\
	memcpy((ptr), in + in_off, (len));	\
	in_off += (len);			\
    } while(0);

fastq *decode_block(unsigned char *in, unsigned int in_size, timings *t) {
    unsigned char *in_end = in + in_size, *comp, *out;
    uint32_t in_off = 0, nr;
    int i, j, err = 0;
    uint32_t u_len, c_len;
    uint8_t c;
    struct timeval tv1, tv2;
    
    GET(&nr, 4);
    fastq *fq = fastq_alloc(nr);

    // ----------
    // Name
    gettimeofday(&tv1, NULL);
    GET(&u_len, 4);
    GET(&c, 1);     // strategy
    GET(&c_len, 4);

    comp = in+in_off;
    in_off += c_len;

    out = (unsigned char *)decode_names(comp, c_len, u_len, c);
    fq->name_buf = (char *)out;
    fq->name_len = u_len;

    gettimeofday(&tv2, NULL);
    t->ncsize += u_len;
    t->nusize += c_len;
    t->ntime += (tv2.tv_sec - tv1.tv_sec) * 1000000;
    t->ntime += tv2.tv_usec - tv1.tv_usec;

    // ----------
    // Lengths
    GET(&c, 1);     // strategy, but also length.  Needed as len?
    if (c > 0) {
	// Fixed length
	uint32_t len;
	int vl;
	in_off += (vl = var_get_u32(in+in_off, in_end, &len));
	err |= vl == 0;

	for (i = 0; i < nr; i++)
	    fq->len[i] = len;
	t->lcsize += nr*4;
	t->lusize += c;
    } else {
	// Variable length
	uint32_t blen;
	GET(&blen, 4); // needed?  Doesn't seem it now
	for (i = 0; i < nr; i++) {
	    int vl;
	    in_off += (vl = var_get_u32(in+in_off, in_end, &fq->len[i]));
	    err |= vl == 0;
	}

	t->lcsize += nr*4;
	t->lusize += blen+5;
    }
    if (err)
	goto err;

    // ----------
    // Seq
    gettimeofday(&tv1, NULL);
    GET(&c, 1);
    GET(&u_len, 4);
    GET(&c_len, 4);
    comp = in+in_off;
    in_off += c_len;

    int slevel = c>>4;
    int both_strands = (c >> 3) & 1;

    if ((c & 7) == 1) {
	out = (uint8_t *)decode_seq(comp, c_len, fq->len, fq->num_records,
				    both_strands, slevel, u_len);
    } else {
	out = rans_uncompress_4x16(comp, c_len, &u_len);
    }

    fq->seq_buf = (char *)out;
    fq->seq_len = u_len;
    for (i = 0; i < nr; i++)
	fq->seq[i] = i ? fq->seq[i-1] + fq->len[i-1] : 0;

    gettimeofday(&tv2, NULL);
    t->stime += (tv2.tv_sec - tv1.tv_sec) * 1000000;
    t->stime += tv2.tv_usec - tv1.tv_usec;
    t->scsize += u_len;
    t->susize += c_len;

    // ----------
    // Qual
    gettimeofday(&tv1, NULL);
    GET(&c, 1);
    GET(&u_len, 4);
    GET(&c_len, 4);
    comp = in+in_off;
    in_off += c_len;

    size_t out_len;
    if (c == 0) {
	// Rans
	out = rans_uncompress_4x16(comp, c_len, &u_len);
	fq->qual_buf = (char *)out;
	fq->qual_len = out_len = u_len;
    } else {
	// FQZComp qual
	fqz_slice s;
	s.num_records = fq->num_records;
	s.len = fq->len;
	s.flags = fq->flag;
	//s.seq = (unsigned char **)fq->seq;
	s.seq = (unsigned char **)malloc(fq->num_records * sizeof(char *));
	for (i = j = 0; i < fq->num_records; j += fq->len[i++])
	    s.seq[i] = (unsigned char *)fq->seq_buf + j;

	// pass lengths as NULL and fix fqz_decompress to cope?
	int *lengths = malloc(nr * sizeof(lengths));
	out = (uint8_t *)fqz_decompress((char *)comp, c_len, &out_len,
					lengths, nr, &s);
	fq->qual_buf = (char *)out;
	fq->qual_len = out_len;
	free(s.seq);
	free(lengths);
    }
    for (i = 0; i < fq->qual_len; i++)
	fq->qual_buf[i] += 33;

    gettimeofday(&tv2, NULL);
    t->qtime += (tv2.tv_sec - tv1.tv_sec) * 1000000;
    t->qtime += tv2.tv_usec - tv1.tv_usec;
    t->qcsize += out_len;
    t->qusize += c_len;

    return fq;

 err:
    free(fq);
    return NULL;
}

typedef struct {
    fqz_gparams *gp;
    opts *arg;
    timings t;
    fastq *fq;
    char *comp;
    uint32_t clen;
    int eof;
} enc_dec_job;

static void *encode_thread(void *arg) {
    enc_dec_job *j = (enc_dec_job *)arg;
    if (j->eof)
	return j;

    j->comp = encode_block(j->gp, j->arg, j->fq, &j->t, &j->clen);
    fastq_free(j->fq);
    return j;
}

#define THREADED

// TODO: use async read and write threads too so main doesn't block on I/O
int encode(FILE *in_fp, FILE *out_fp, fqz_gparams *gp, opts *arg,
	   timings *t) {
    int rans_methods = (1<<RANS0) | (1<<RANS1) | (1<<RANS129) | (1<<RANS193);
    // Name
    method_avail[SEC_NAME] = 1<<LZP3;
    if (arg->nstrat == 1) {
	method_avail[SEC_NAME] |= 1<<(TOK3_3_LZP + arg->nlevel/2-1);
	method_avail[SEC_NAME] |= 1<<(TOK3_3 + arg->nlevel/2-1);
    }
    if (arg->nstrat == 2)
	method_avail[SEC_NAME] |= 1<<(TOK3_3_LZP + arg->nlevel/2-1);
    
    // Seq
    method_avail[SEC_SEQ] = rans_methods;
    if (arg->sstrat == 1) {
	if (arg->slevel >= 10)
	    method_avail[SEC_SEQ] |= 1<<SEQ10;
	if (arg->slevel >= 12) {
	    if (arg->both_strands | arg->slevel > 12)
		method_avail[SEC_SEQ] |= 1<<SEQ12B;
	    else
		method_avail[SEC_SEQ] |= 1<<SEQ12;
	}
	if (arg->slevel >= 13)
	    method_avail[SEC_SEQ] |= 1<<SEQ13B;
	if (arg->slevel >= 14)
	    method_avail[SEC_SEQ] |= 1<<SEQ14B;
    }

    // Qual
    method_avail[SEC_QUAL] = rans_methods;
    if (arg->qstrat == 1) {
	int m;
	for (m = 0; m <= arg->qlevel && m <= 3; m++)
	    method_avail[SEC_QUAL] |= 1<<(FQZ0+m);
    }

#ifdef THREADED
    int n = arg->nthread, end = 0;
    hts_tpool *p = hts_tpool_init(n);
    hts_tpool_process *q = hts_tpool_process_init(p, n*2, 0);
    hts_tpool_result *r;
    enc_dec_job *j;
#endif

    char *in = malloc(arg->blk_size);
    int in_rem = 0;

    for(;;) {
	int nbytes = fread(in+in_rem, 1, arg->blk_size - in_rem, in_fp);
	if (nbytes < 0)
	    return -1;
	if (nbytes == 0)
	    break;
	nbytes += in_rem;

	fastq *fq = load_seqs(in, nbytes, &in_rem);
	memmove(in, in+in_rem, nbytes - in_rem);
	in_rem = nbytes - in_rem;

	if (!fq)
	    return -1;

	if (!fq->num_records) {
	    fastq_free(fq);
	    break;
	}

#ifdef THREADED
	// Dispatch a job
	j = calloc(1, sizeof(*j));
	j->gp = gp;
	j->arg = arg;
	memset(&j->t, 0, sizeof(j->t));
	j->fq = fq;
	j->eof = 0;

	if (hts_tpool_dispatch(p, q, encode_thread, j) != 0)
	    goto err;

	// Check for a result
	if ((r = hts_tpool_next_result(q))) {
	    j = hts_tpool_result_data(r);
	    if (j->eof) {
		end = 1;
	    } else {
		append_timings(t, &j->t, arg->verbose);
		fwrite(&j->clen, 1, 4, out_fp);
		fwrite(j->comp, 1, j->clen, out_fp);
		free(j->comp);
	    }
	    hts_tpool_delete_result(r, 1);
	}
#else
	uint32_t clen;
	t->nblock++;
	char *out = encode_block(gp, arg, fq, t, &clen);
	fwrite(&clen, 1, 4, out_fp);
	fwrite(out, 1, clen, out_fp);
	free(out);

	fastq_free(fq);
#endif
    }

#ifdef THREADED
    j = malloc(sizeof(*j));
    j->eof = 1;
    if (hts_tpool_dispatch(p, q, encode_thread, j) != 0)
	goto err;

    // End of input, so work through remaining results
    while (!end && (r = hts_tpool_next_result_wait(q))) {
	enc_dec_job *j = hts_tpool_result_data(r);
	if (j->eof) {
	    end = 1;
	} else {
	    append_timings(t, &j->t, arg->verbose);
	    fwrite(&j->clen, 1, 4, out_fp);
	    fwrite(j->comp, 1, j->clen, out_fp);
	    free(j->comp);
	}
	hts_tpool_delete_result(r, 1);
    }
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);
#endif

    return 0;

 err:
    return -1;
}

int output_fastq(FILE *out_fp, fastq *fq) {
    char *np = fq->name_buf;
    char *sp = fq->seq_buf;
    char *qp = fq->qual_buf;

#if 1
    // A bit faster sometimes.
    int len = fq->name_len + fq->seq_len + fq->qual_len + fq->num_records*5;
    char *buf = malloc(len), *cp = buf;

    for (int i = 0; i < fq->num_records; i++) {
	*cp++ = '@';
	while ((*cp++ = *np++))
	    ;
	*--cp = '\n'; cp++;
	memcpy(cp, sp, fq->len[i]);
	cp += fq->len[i];
	sp += fq->len[i];
	*cp++ = '\n';
	*cp++ = '+';
	*cp++ = '\n';
	memcpy(cp, qp, fq->len[i]);
	cp += fq->len[i];
	qp += fq->len[i];
	*cp++ = '\n';
    }
    fwrite(buf, 1, cp-buf, out_fp);
    free(buf);
#else
    for (int i = 0; i < fq->num_records; i++) {
	fprintf(out_fp, "@%s\n%.*s\n+\n%.*s\n",
		np, fq->len[i], sp, fq->len[i], qp);
	np += strlen(np)+1;
	sp += fq->len[i];
	qp += fq->len[i];
    }
#endif

    return 0;
}

static void *decode_thread(void *arg) {
    enc_dec_job *j = (enc_dec_job *)arg;
    if (j->eof)
	return j;

    j->fq = decode_block(j->comp, j->clen, &j->t);
    free(j->comp);
    return j;
}

int decode(FILE *in_fp, FILE *out_fp, opts *arg, timings *t) {
#ifdef THREADED
    int n = arg->nthread, end = 0;
    hts_tpool *p = hts_tpool_init(n);
    hts_tpool_process *q = hts_tpool_process_init(p, n*2, 0);
    hts_tpool_result *r;
    enc_dec_job *j;
#endif

    for (;;) {
	// Load next compressed block
	int i, c_len;
	unsigned char *comp;

	// Total compressed size
	if (fread(&c_len, 1, 4, in_fp) != 4)
	    break;

	comp = malloc(c_len);
	if (fread(comp, 1, c_len, in_fp) != c_len)
	    return -1;

#ifdef THREADED
	// Dispatch a job
	j = calloc(1, sizeof(*j));
	memset(&j->t, 0, sizeof(j->t));
	j->comp = comp;
	j->clen = c_len;
	j->fq = NULL;
	j->eof = 0;

	if (hts_tpool_dispatch(p, q, decode_thread, j) != 0)
	    goto err;

	// Check for a result
	if ((r = hts_tpool_next_result(q))) {
	    j = hts_tpool_result_data(r);
	    if (j->eof) {
		end = 1;
	    } else {
		append_timings(t, &j->t, arg->verbose);
		output_fastq(out_fp, j->fq);
		free(j->fq);
	    }
	    hts_tpool_delete_result(r, 1);
	}
#else
	t->nblock++;
	fastq *fq = decode_block(comp, c_len, t);

	// ----------
	// Convert back to fastq
	char *np = fq->name_buf;
	char *sp = fq->seq_buf;
	char *qp = fq->qual_buf;

	for (i = 0; i < fq->num_records; i++) {
	    fprintf(out_fp, "@%s\n%.*s\n+\n%.*s\n",
		    np, fq->len[i], sp, fq->len[i], qp);
	    np += strlen(np)+1;
	    sp += fq->len[i];
	    qp += fq->len[i];
	}

	fastq_free(fq);
	free(comp);
#endif
    }

#ifdef THREADED
    j = malloc(sizeof(*j));
    j->eof = 1;
    if (hts_tpool_dispatch(p, q, decode_thread, j) != 0)
	goto err;

    // End of input, so work through remaining results
    while (!end && (r = hts_tpool_next_result_wait(q))) {
	enc_dec_job *j = hts_tpool_result_data(r);
	if (j->eof) {
	    end = 1;
	} else {
	    append_timings(t, &j->t, arg->verbose);
	    output_fastq(out_fp, j->fq);
	    free(j->fq);
	}
	hts_tpool_delete_result(r, 1);
    }
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);
#endif

    return 0;

 err:
    // FIXME: tidy up
    return -1;
}

void usage(FILE *fp) {
    fprintf(fp, "Usage: fqzcomp5 [options]    [input.fastq [output.fqz5]]\n");
    fprintf(fp, "Usage: fqzcomp5 [options] -d [input.fqz5  [output.fastq]]\n");
    fprintf(fp, "\nOptions:\n");
    fprintf(fp, "    -d            Decompress\n");
    fprintf(fp, "    -t INT        Number of threads.  Defaults to 4\n");
    fprintf(fp, "    -b SIZE       Specify block size. May use K, M and G sufixes\n");
    fprintf(fp, "    -v            Increase verbostity\n");
    fprintf(fp, "    -V            Silent mode\n");
    fprintf(fp, "\n");
    fprintf(fp, "    -n INT        Name encoding method (0=rANS, 1=tok3, 2=tok3+LZP)\n");
    fprintf(fp, "    -N INT        Name encoding strategy.\n");
    fprintf(fp, "    -s INT        Sequence encoding method (0=rANS, 1=fqz)\n");
    fprintf(fp, "    -S INT        Sequence encoding strategy (context size)\n");
    fprintf(fp, "    -B            Update sequence context on both strands\n");
    fprintf(fp, "    -q INT        Quality encoding method (0=rANS, 1=fqz)\n");
    fprintf(fp, "    -Q INT        Quality encoding strategy (0 to 3)\n");
    fprintf(fp, "\n");
    fprintf(fp, "Compression levels:\n");
    fprintf(fp, "    -1            Equivalent to -n0 -s0 -q0 -b10M\n");
    fprintf(fp, "    -3            Equivalent to -n1 -s0 -q0 -b100M\n");
    fprintf(fp, "    -5            Equivalent to -n2 -s1 -q1 -b100M\n");
    fprintf(fp, "    -7            Equivalent to -n2 -s1 -q1 -b500M -B -S14\n");
    fprintf(fp, "    -9            Equivalent to -n2 -s1 -q1 -b1GM  -B -S15\n");
}

int main(int argc, char **argv) {
    int decomp = 0;
    fqz_gparams *gp = NULL, gp_local;
    FILE *in_fp, *out_fp;
    timings t = {0};

    opts arg = {
	.qstrat = 1, // 0=rans, 1=fqz
	.qlevel = 0,
	.sstrat = 1, // 0=rans, 1=fqz
	.slevel = 12,// seq context = 4^12 
	.nstrat = 2, // (0=rans), 1=tok3, 2=tok3 + comments
	.nlevel = 5,
	.both_strands =0, // adjusts seq strat 1.
	.verbose = 0,
	.blk_size = BLK_SIZE,
	.nthread = 4,
    };

#ifdef _WIN32
        _setmode(_fileno(stdin),  _O_BINARY);
        _setmode(_fileno(stdout), _O_BINARY);
#endif

    extern char *optarg;
    extern int optind;
    int opt;

    while ((opt = getopt(argc, argv, "dq:Q:b:x:Bs:S:vn:N:Vt:h13579")) != -1) {
	switch (opt) {
	case 't':
	    arg.nthread = atoi(optarg);
	    if (arg.nthread < 1)
		arg.nthread = 1;
	    break;

	case 'v':
	    arg.verbose++;
	    break;

	case 'V':
	    arg.verbose = -1;
	    break;

	case 'd':
	    decomp = 1;
	    break;

	case 'B':
	    arg.both_strands=1;
	    break;

	case 's':
	    arg.sstrat = atoi(optarg);
	    break;
	case 'S':
	    arg.slevel = atoi(optarg);
	    if (arg.slevel < 0)
		arg.slevel = 0;
	    if (arg.slevel > 16)
		arg.slevel = 16;
	    break;

	case 'n':
	    arg.nstrat = atoi(optarg);
	    break;
	case 'N':
	    arg.nlevel = atoi(optarg);
	    if (arg.nlevel < 0)
		arg.nlevel = 0;
	    if (arg.nlevel > 19)
		arg.nlevel = 19;
	    break;

	case 'q':
	    arg.qstrat = atoi(optarg);
	    break;
	case 'Q':
	    arg.qlevel = atoi(optarg);
	    break;

	case 'b': {
	    char *endp;
	    arg.blk_size = strtol(optarg, &endp, 0);
	    if (*endp == 'k' || *endp == 'K')
		arg.blk_size *= 1000;
	    else if (*endp == 'm' || *endp == 'M')
		arg.blk_size *= 1000000;
	    else if (*endp == 'g' || *endp == 'G')
		arg.blk_size *= 1000000000;
	    if (arg.blk_size < 1000000)
		arg.blk_size = 1000000;
	    if (arg.blk_size > 2000000000)
		arg.blk_size = 2000000000;
	    break;
	}

	case '1':
	    arg.nstrat = 0;
	    arg.sstrat = 0;
	    arg.qstrat = 0;
	    arg.blk_size = 10e6;
	    break;

	case '3':
	    arg.nstrat = 1;
	    arg.nlevel = 3;
	    arg.sstrat = 0;
	    arg.qstrat = 0;
	    arg.blk_size = 100e6;
	    break;

	case '5':
	    arg.nstrat = 2;
	    arg.nlevel = 5;
	    arg.sstrat = 1;
	    arg.qstrat = 1;
	    arg.blk_size = 100e6;
	    arg.qlevel = 1;
	    break;

	case '7':
	    // TODO: also add format detection, so qlevel is adjusted
	    // per file.  Or auto-sensing and learning strategy (ala CRAM)
	    arg.nstrat = 2;
	    arg.nlevel = 7;
	    arg.sstrat = 1;
	    arg.both_strands = 1;
	    arg.slevel = 14;
	    arg.qstrat = 1;
	    arg.qlevel = 3;
	    break;

	case '9':
	    arg.nstrat = 2;
	    arg.nlevel = 9;
	    arg.sstrat = 1;
	    arg.both_strands = 1;
	    arg.slevel = 15;
	    arg.qstrat = 1;
	    arg.qlevel = 3;
	    arg.blk_size = 1e9;
	    break;

#if 0
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
#endif

	case 'h':
	    usage(stdout);
	    return 0;

	default:
	    usage(stderr);
	    return 1;
	}
    }

    if (optind == argc && isatty(0)) {
	usage(stdout);
	return 0;
    }

    in_fp = optind < argc ? fopen(argv[optind], "r") : stdin;
    out_fp = ++optind < argc ? fopen(argv[optind], "wb") : stdout;

    // FIXME: use variable sized integers

    // Block based, for arbitrary sizes of input
    if (decomp) {
	if (decode(in_fp, out_fp, &arg, &t) < 0)
	    exit(1);
    } else {
	if (encode(in_fp, out_fp, gp, &arg, &t) < 0)
	    exit(1);
    }

    if (arg.verbose >= 0) {
	fprintf(stderr, "All %ld blocks combined:\n", t.nblock);
	fprintf(stderr, "Names    %10ld to %10ld in %.2f sec\n",
		t.nusize, t.ncsize, t.ntime/1e6);
	fprintf(stderr, "Lengths  %10ld to %10ld\n",
		t.lusize, t.lcsize);
	fprintf(stderr, "Seqs     %10ld to %10ld in %.2f sec\n", 
		t.susize, t.scsize, t.stime/1e6);
	fprintf(stderr, "Qual     %10ld to %10ld in %.2f sec\n", 
		t.qusize, t.qcsize, t.qtime/1e6);
    }

    fclose(in_fp);
    fclose(out_fp);

    return 0;
}
