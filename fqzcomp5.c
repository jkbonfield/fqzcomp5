// TODO
// - Tokenise name with different alphabets.
//   See ~/scratch/data/enano_vir/ERR2708432.fastq

// - Split aux tags into own data series using CRAM's TL + per tag.

// - Seq encoding using STR + copy-number?
//   May work well for ONT/PacBio where run len varies.

// - Also permit old name encoding (fqzcomp -n1).
//   Better on above data. Why?
//   (fqz -n1: 912989, -n2: 1340223, paq8: 891928

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
#include <sys/time.h>

#include "htscodecs/varint.h"
#include "htscodecs/fqzcomp_qual.h"
#include "htscodecs/tokenise_name3.h"
#include "htscodecs/rANS_static4x16.h"
#include "htscodecs/varint.h"
#include "lzp16e.h"

#define BLK_SIZE 512*1000000

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
		 int *len, int nrecords, int both_strands, int ctx_size,
		 /*unsigned char *out,*/ unsigned int *out_size) {
    char *out = malloc(in_size + 100);
    if (!out)
	return NULL;

    const int msize = 1<<(2*ctx_size);
    const int mask = msize-1;
    int i;

    SMALL_MODEL(4,_) *seq_model = malloc(msize * sizeof(*seq_model));
//    SMALL_MODEL(4,_) *hom_len = malloc(msize * sizeof(*hom_len));
    for (i = 0; i < msize; i++) {
	SMALL_MODEL(4,_init)(&seq_model[i]);
//	SMALL_MODEL(4,_init)(&hom_len[i]);
    }

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
	int rep_len = 0;
	int last_base = 99;
	int last_base2 = 98;
	switch (state) {
	case uc_ACGT:
	case lc_ACGT:
	    for (j = 0; j < run; j++) {
#if 0
		// Encode GAAAT as G0 A3 T0
		int k;
		for (k = j+1; k < run; k++) {
		    if (in[i+k] != in[i+j])
			break;
		}
		int r2 = k-j-1;
		//int rctx = ((last & 0x3ff)<<4);
		int rctx = last & (mask>>2);
		int n = 0;
		do {
		    SMALL_MODEL(4,_encodeSymbol)(&hom_len[rctx & 0xfff], &rc,
						 MIN(3,r2));
		    rctx = (rctx & (mask>>2)) | (n<<(2*ctx_size-2));
		    r2 -= 3;
		    n+=(n<3);
		} while(r2 >= 0);
 		j = k-1;
#endif

// Use previous repeat length as part of seq context.
// Also poor (but not as bad as above)
//		rep_len = (last_base == last_base2)
//		    ? rep_len+(rep_len<3)
//		    : 0;
//		last_base2 = last_base;
//		last_base = in[i+j];
//		unsigned char b = L[in[i+j]] & 3;
//		int ctx = (last&(mask>>2)) | (rep_len<<(2*ctx_size-2));
//		SMALL_MODEL(4, _encodeSymbol)(&seq_model[ctx], &rc, b);

		unsigned char b = L[in[i+j]] & 3;
		SMALL_MODEL(4, _encodeSymbol)(&seq_model[last], &rc, b);

		last = ((last<<2) + b) & mask;
		_mm_prefetch((const char *)&seq_model[(last<<4)&mask],
			     _MM_HINT_T0);

		// 0.7% and 3.2% smaller for _.FQ and _.fq respectively (at ctx_size 12),
		// but 45% more CPU for seq encoding.
		if (both_strands) {
		    int b2 = last2 & 3;
		    last2 = last2/4 + ((3-b) << (2*ctx_size-2));
		    // Can't predict prefetch for other end.
		    //_mm_prefetch((const char *)&seq_model[last2], _MM_HINT_T0);
		    SMALL_MODEL(4, _updateSymbol)(&seq_model[last2], b2);
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
		 int *len, int nrecords, int both_strands, int ctx_size,
		 /*unsigned char *out,*/ unsigned int out_size) {
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
    RC_SetInput(&rc, in, in+in_size);
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

typedef struct {
    int nstrat, sstrat, qstrat;
    int nlevel, slevel, qlevel;
    int verbose;
    int both_strands;
    int blk_size;
} opts;

typedef struct {
    int64_t nusize, ncsize, ntime;
    int64_t susize, scsize, stime;
    int64_t qusize, qcsize, qtime;
    int64_t lusize, lcsize, ltime;
    int64_t osize;
} timings;

int encode(FILE *in_fp, FILE *out_fp, fqz_gparams *gp, opts *arg, timings *t) {
    unsigned char *out;
    size_t out_len;
    struct timeval tv1, tv2;
    int c, err = 0;

    for(;;) {
	fastq *fq = load_seqs(in_fp, arg->blk_size);
	if (!fq)
	    return -1;

	if (arg->verbose)
	    fprintf(stderr, "Loaded nrec = %d\n", fq->num_records);
	if (!fq->num_records) {
	    fastq_free(fq);
	    break;
	}
	//fastq_dump(fq);
	    
	fwrite(&fq->num_records, 1, 4, out_fp);

	//----------
	// Names: tok3
	// Strat 0 = LZP + rANS
	// Strat 1 = Tok3
	// Strat 2 = Name(tok3)+Flag(RC)+Comment(LZP+rANS)
	gettimeofday(&tv1, NULL);
	int clen;

	fwrite(&fq->name_len, 1, 4, out_fp);
	if (arg->nstrat == 0) {
	    // TODO: work out a better maximum bound
	    char *lzp_out = malloc(fq->name_len*2);
	    clen = lzp(fq->name_buf, fq->name_len, lzp_out);
	    out = rans_compress_4x16(lzp_out, clen, &clen, 5);
	    free(lzp_out);
	    putc(0, out_fp); // name method
	    fwrite(&clen, 1, 4, out_fp);
	    fwrite(out, 1, clen, out_fp);
	    t->osize += 9;
	    free(out);

	} else if (arg->nstrat == 1) {
	    out = tok3_encode_names(fq->name_buf, fq->name_len, arg->nlevel,
		                    0, &clen, NULL);
	    putc(1, out_fp); // name method
	    fwrite(&clen, 1, 4, out_fp);
	    fwrite(out, 1, clen, out_fp);
	    t->osize += 9;
	    free(out);

	} else {
	    char *n1 = malloc(fq->name_len);
	    char *n2 = malloc(fq->name_len);
	    char *flag = malloc(fq->name_len/2);  //Worst case\n 
	    char *cp1 = n1, *cp2 = n2;
	    int i = 0, nr = 0;
	    // Flag bit 0: has "/NUM"
	    // Flag bit 1: /1 vs /2
	    // Flag bit 2: has a comment
	    // Flag bit 3: space vs tab before comment
	    while (i < fq->name_len) {
		int j, k, f = 0;
		int w1end = 0;
		int w2start = 0;
		int w2end = 0;
		for (j = i; j < fq->name_len; j++) {
                   if (fq->name_buf[j] == '\0') {
		       w2end = j;
		       break;
	           }
		   if (!w2start && (fq->name_buf[j] == ' ' ||
                                    fq->name_buf[j] == '\t')) {
		       w2end = j;
		       w2start = j+1;
		       f |= 4; // FLAG: has comment
	           }
	        }

		if (!w2end)
		    w2end = fq->name_len;

		if (w2start)
		    // FLAG: space vs tab
		    f |= fq->name_buf[w2start-1] == ' ' ? 0 : 8;

		if (w2end>1 && fq->name_buf[w2end-2] == '/') {
		    // FLAG /1 or /2
		    if (fq->name_buf[w2end-1] == '1')
			f |= 1, w2end -= 2;
		    else if (fq->name_buf[w2end-1] == '2')
			f |= 3, w2end -= 2;
	        }

		flag[nr++] = f;
		memcpy(cp1, &fq->name_buf[i], w2end-i);
		cp1[w2end-i]=0;
		cp1 += w2end-i+1;

		if (w2start) {
		    memcpy(cp2, &fq->name_buf[w2start], w2end-w2start);
		    cp2[w2end-w2start] = 0;
		    cp2 += w2end-w2start+1;
		}

		i = j+1;
	    }

	    int clen1, clen2 = 0, clenf;
	    out = tok3_encode_names(n1, cp1-n1, arg->nlevel, 0, &clen1, NULL);
	    char *outf = rans_compress_4x16(flag, nr, &clenf, 129);
	    char *out2 = NULL;
	    if (cp2 != n2) {
		char *lzp_out = malloc((cp2-n2)*2);
		clen2 = lzp(n2, cp2-n2, lzp_out);
		out2 = rans_compress_4x16(lzp_out, clen2, &clen2, 5);
		free(lzp_out);
	    }

	    clen = clen1 + clenf + clen2 + 8;

	    putc(2, out_fp); // name method
	    fwrite(&clen,  1, 4, out_fp);
	    fwrite(&clen1, 1, 4, out_fp);
	    fwrite(&clenf, 1, 4, out_fp);

	    fwrite(out, 1, clen1, out_fp);
	    free(out);

	    fwrite(outf, 1, clenf, out_fp);
	    free(outf);

	    if (out2) {
		fwrite(out2, 1, clen2, out_fp);
		free(out2);
	    }
	    t->osize += 12; // should add to clen instead?
	    free(n1);
	    free(n2);
	    free(flag);
	}
	if (arg->verbose)
	    fprintf(stderr, "Names: %10d to %10d\n", fq->name_len, clen);
	t->nusize += fq->name_len;
	t->ncsize += clen;
	gettimeofday(&tv2, NULL);
	t->ntime += (tv2.tv_sec - tv1.tv_sec) * 1000000;
	t->ntime += tv2.tv_usec - tv1.tv_usec;

	//----------
	// Read lengths
	if (fq->fixed_len) {
	    // Fixed length, with next byte holding the size of length
	    char buf[5], nb = 1;
	    nb += var_put_u32(buf+1, NULL, fq->fixed_len);
	    buf[0] = nb-1;
	    t->osize += 0;
	    t->lusize += 4*fq->num_records;
	    t->lcsize += nb;
	    fwrite(buf, 1, nb, out_fp);
	    if (arg->verbose)
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
	    t->osize += 0;
	    t->lusize += fq->num_records*4;
	    t->lcsize += nb+5;
	    fwrite(buf, 1, nb, out_fp);
	    free(buf);
	    if (arg->verbose)
		fprintf(stderr, "Len:   %10d to %10d\n",
			4*fq->num_records, nb+5);
	}

	//----------
	// Seq: rans or statistical modelling
	gettimeofday(&tv1, NULL);
	if (arg->sstrat == 1) {
	    out = encode_seq(fq->seq_buf, fq->seq_len,
			     fq->len, fq->num_records,
			     arg->both_strands, arg->slevel, &clen);
	    if (!out) {
		fprintf(stderr, "ERR: failed to encode sequence\n");
		return -1;
	    }
	    putc((arg->slevel<<4) | (arg->both_strands<<3) | 1, out_fp);
	    fwrite(&fq->seq_len, 1, 4, out_fp);
	    fwrite(&clen, 1, 4, out_fp);
	    fwrite(out, 1, clen, out_fp);
	    t->osize += 9;
	} else {
	    // FIXME: test 197, 65 and 1
	    // FIXME: also try fqz for encoding, qmap and qshift=2.
	    // But need extension to permit reverse complement.
	    // Also not quite the same in terms of model updates?
	    out = rans_compress_4x16(fq->seq_buf, fq->seq_len, &clen, 197);
	    putc(0, out_fp); // seq method
	    fwrite(&fq->seq_len, 1, 4, out_fp);
	    fwrite(&clen, 1, 4, out_fp);
	    fwrite(out, 1, clen, out_fp);
	    t->osize += 9;
	}
	free(out);

	if (arg->verbose)
	    fprintf(stderr, "Seq:   %10d to %10d\n", fq->seq_len, clen);
	t->susize += fq->seq_len;
	t->scsize += clen;
	gettimeofday(&tv2, NULL);
	t->stime += (tv2.tv_sec - tv1.tv_sec) * 1000000;
	t->stime += tv2.tv_usec - tv1.tv_usec;

	//----------
	// Qual: rans or fqz
	// Convert fastq struct to fqz_slice for context
	gettimeofday(&tv1, NULL);
	if (arg->qstrat == 0) {
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
	    out = fqz_compress(4, s, fq->qual_buf, fq->qual_len,
			       &out_len, arg->qlevel, gp);
	    free(s->seq);
	    free(s);
	}
	fwrite(&fq->qual_len, 1, 4, out_fp);
	fwrite(&out_len, 1, 4, out_fp);
	fwrite(out, 1, out_len, out_fp);
	t->osize += 9;
	free(out);

	if (arg->verbose)
	    fprintf(stderr, "Qual:  %10d to %10d\n", fq->qual_len, (int)out_len);
	t->qusize += fq->qual_len;
	t->qcsize += out_len;
	gettimeofday(&tv2, NULL);
	t->qtime += (tv2.tv_sec - tv1.tv_sec) * 1000000;
	t->qtime += tv2.tv_usec - tv1.tv_usec;

	fastq_free(fq);
    }

    return 0;
}

int decode(FILE *in_fp, FILE *out_fp, opts *arg, timings *t) {
    unsigned char *out;
    size_t out_len;
    struct timeval tv1, tv2;
    int c, err = 0;

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
	if ((c = getc(in_fp)) == EOF) break;
	if (fread(&c_len, 1, 4, in_fp) != 4)
	    break;
	comp = malloc(c_len);
	if (fread(comp, 1, c_len, in_fp) != c_len)
	    break;

	gettimeofday(&tv1, NULL);
	if (c == 0) {
	    int ru_len;
	    char *rout = rans_uncompress_4x16(comp, c_len, &ru_len);
	    out = malloc(u_len);
	    u_len = unlzp(rout, ru_len, out);
	    free(rout);
	} else if (c == 1) {
	    out = tok3_decode_names(comp, c_len, &u_len);
	} else {
	    uint32_t clen1 = *(uint32_t *)comp;
	    uint32_t clenf = *(uint32_t *)(comp+4);
	    uint32_t clen2 = c_len - clen1 - clenf - 8;

	    // Uncompress 3 separate components
	    int u_len1, u_lenf, u_len2, ru_len;
	    char *out1 = tok3_decode_names(comp+8, clen1, &u_len1);
	    char *outf = rans_uncompress_4x16(comp+8+clen1, clen1, &u_lenf);
	    char *out2 = NULL;
	    if (clen2) {
		int rulen;
		char *rout = rans_uncompress_4x16(comp+8+clen1+clenf, clen2,
		                                  &rulen);
		out2 = malloc(u_len);
		u_len2 = unlzp(rout, rulen, out2);
		free(rout);
	    }

	    // Stitch together ID + flag + comment
	    char *cp1 = out1, *cp1_end = out1+u_len1;
	    char *cpf = outf, *cpf_end = outf+u_lenf;
	    char *cp2 = out2, *cp2_end = out2 + u_len2;
	    out = malloc(u_len);
	    char *cp  = out,  *cp_end = out + u_len;
	    char *last_cp = NULL;
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
		
		if (cp2) {
		    while (cp2 < cp2_end && cp < cp_end && *cp2)
			*cp++ = *cp2++;
		    cp2++;
	        }

		if (cp == last_cp)
		    // ran out of data early; avoids looping forever
		    break;

		*cp++ = 0;
		last_cp = cp;
	    }

	    free(out1);
	    free(outf);
	    free(out2);
	}
	fq->name_buf = out;
	fq->name_len = u_len;
	free(comp);
	gettimeofday(&tv2, NULL);
	t->ncsize += u_len;
	t->nusize += c_len;
	t->ntime += (tv2.tv_sec - tv1.tv_sec) * 1000000;
	t->ntime += tv2.tv_usec - tv1.tv_usec;

	// ----------
	// Lengths
	if ((c = getc(in_fp)) == EOF) break;
	if (c > 0) {
	    // Fixed length
	    char buf[10];
	    err |= fread(buf, 1, c, in_fp) != c;
	    uint32_t len;
	    err |= var_get_u32(buf, buf+c, &len) == 0;
	    for (i = 0; i < nr; i++)
		fq->len[i] = len;
	    t->lcsize += nr*4;
	    t->lusize += c+4;
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
	    t->lcsize += nr*4;
	    t->lusize += blen+4;
	}
	if (err)
	    break;

	// ----------
	// Seq
	if ((c = getc(in_fp)) == EOF) break;
	arg->slevel = c>>4;
	arg->both_strands = (c >> 3) & 1;
	if (fread(&u_len, 1, 4, in_fp) != 4)
	    break;
	if (fread(&c_len, 1, 4, in_fp) != 4)
	    break;
	comp = malloc(c_len);
	if (fread(comp, 1, c_len, in_fp) != c_len)
	    break;
	gettimeofday(&tv1, NULL);
	if ((c & 7) == 1) {
	    out = decode_seq(comp, c_len, fq->len, fq->num_records,
			     arg->both_strands, arg->slevel, u_len);
	} else {
	    out = rans_uncompress_4x16(comp, c_len, &u_len);
	}
	free(comp);
	fq->seq_buf = out;
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
	int mode = getc(in_fp);
	if (fread(&u_len, 1, 4, in_fp) != 4)
	    break;
	if (fread(&c_len, 1, 4, in_fp) != 4)
	    break;
	comp = malloc(c_len);
	if (fread(comp, 1, c_len, in_fp) != c_len)
	    break;

	gettimeofday(&tv1, NULL);
	if (mode == 0) {
	    // Rans
	    out = rans_uncompress_4x16(comp, c_len, &u_len);
	    fq->qual_buf = out;
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
	gettimeofday(&tv2, NULL);
	t->qtime += (tv2.tv_sec - tv1.tv_sec) * 1000000;
	t->qtime += tv2.tv_usec - tv1.tv_usec;
	t->qcsize += out_len;
	t->qusize += c_len;

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

    return 0;
}

int main(int argc, char **argv) {
    unsigned char *in, *out, *seq;
    size_t in_len, out_len, seq_len;
    int decomp = 0, vers = 4;  // CRAM version 4.0 (4) or 3.1 (3)
    fqz_gparams *gp = NULL, gp_local;
    FILE *in_fp, *out_fp;
    timings t = {0};

    opts arg = {
	.qstrat = 1, // 0=rans, 1=fqz
	.qlevel = 0,
	.sstrat = 1, // 0=rans, 1=fqz
	.slevel = 12,
	.nstrat = 1, // (0=rans), 1=tok3
	.nlevel = 7,
	.both_strands =0,
	.verbose = 0,
	.blk_size = BLK_SIZE,
    };

#ifdef _WIN32
        _setmode(_fileno(stdin),  _O_BINARY);
        _setmode(_fileno(stdout), _O_BINARY);
#endif

    extern char *optarg;
    extern int optind;
    int opt;

    while ((opt = getopt(argc, argv, "dq:Q:b:x:Bs:S:vn:N:V")) != -1) {
	switch (opt) {
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
	    if (arg.blk_size < 100000)
		arg.blk_size = 100000;
	    break;
	}

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
	}
    }

    in_fp = optind < argc ? fopen(argv[optind], "r") : stdin;
    out_fp = ++optind < argc ? fopen(argv[optind], "wb") : stdout;

    // FIXME: use variable sized integers

    // Block based, for arbitrary sizes of input
    int64_t ntime = 0, stime = 0, qtime = 0;
    size_t nusize = 0, ncsize = 0;
    size_t lusize = 0, lcsize = 0;
    size_t susize = 0, scsize = 0;
    size_t qusize = 0, qcsize = 0;
    size_t osize = 0;
    struct timeval tv1, tv2;
    if (decomp) {
	if (decode(in_fp, out_fp, &arg, &t) < 0)
	    exit(1);
    } else {
	if (encode(in_fp, out_fp, gp, &arg, &t) < 0)
	    exit(1);
    }

    if (arg.verbose >= 0) {
	fprintf(stderr, "Name:    %10ld to %10ld in %.2f sec\n",
		t.nusize, t.ncsize, t.ntime/1e6);
	fprintf(stderr, "Length:  %10ld to %10ld\n",
		t.lusize, t.lcsize);
	fprintf(stderr, "Seq:     %10ld to %10ld in %.2f sec\n", 
		t.susize, t.scsize, t.stime/1e6);
	fprintf(stderr, "Qual:    %10ld to %10ld in %.2f sec\n", 
		t.qusize, t.qcsize, t.qtime/1e6);
	fprintf(stderr, "Other:   %10ld\n",
		t.osize);
    }

    fclose(in_fp);
    fclose(out_fp);

    return 0;
}
