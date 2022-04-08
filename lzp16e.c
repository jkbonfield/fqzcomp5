/*
 * Copyright (c) 2021-2022 Genome Research Ltd.
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
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

// HASH_LEN dictates memory usage and speed.
// It needs 4*2^HASH_LEN bytes for the table.
// Changing this requires rebuilding the decoder too, so it
// needs to be part of the format or cast in stone.
#ifndef HASH_LEN
#    define HASH_LEN 16
#endif

// MIN_LEN dictates when we bother to encode an LZP match.
// Changing this does not affect the format and the decoder can always
// read the data regardless of MIN_LEN.
#ifndef MIN_LEN
#    define MIN_LEN 3
#endif

#define MATCH_CHAR 233
//#define MATCH_CHAR 'Y' // deliberately something in our test data

// Compare &in[i] to m and return match length.
static int match_len(unsigned char *in, int i, int in_len, unsigned char *m) {
    int ml = 0;

    // Safe without unaligned accesses, but fast enough
    in_len -= i;
    in     += i;

    // Rapid discard if it doesn't pass the min len check
    if (in_len < MIN_LEN || memcmp(in, m, MIN_LEN))
	return 0;

#if 1
    // You'd think this would speed it up, but not reliably so!
    if (in_len > MIN_LEN) {
	m  += MIN_LEN;
	in += MIN_LEN;
	ml += MIN_LEN;
    }

    while (ml < in_len && *in == *m) {
	in++;
	m++;
	ml++;
    }
#else
    ml = MIN_LEN;
    while (ml < in_len-16 && memcmp(in+ml, m+ml, 16) == 0)
	ml+=16;
    if (ml < in_len-16)
	while (in[ml] == m[ml])
	    ml++;
    else
	while (ml < in_len && in[ml] == m[ml])
	    ml++;
#endif
    
    return ml;
}

//#define update_hash(h,c) ((((h)<<8) + (c)) & ((1<<HASH_LEN)-1))
//#define update_hash(h,c) ((((h)<<5) ^ (c)) & ((1<<HASH_LEN)-1))
//#define update_hash(h,c) ((((h)<<4) ^ (c)) & ((1<<HASH_LEN)-1))
//#define update_hash(h,c) ((((h*147483475)<<4) ^ c) & ((1<<HASH_LEN)-1))
//#define update_hash(h,c) ((((h)<<4) + (c)*147483475) & ((1<<HASH_LEN)-1))
#define update_hash(h,c) (((((h*0x8ca6b53)<<4)+(h<<5)*17) ^ c) & ((1<<HASH_LEN)-1))

/*
With lzp we have:
in[]  - our input buffer, and a position within it in[i].
h     - a hash key based on previous data e.g. in[i-1].
ht[]  - hash table indexed by h pointing to next bytes of in[],
        for the most recent location matched in 'in'.

Thus ht[h] predicts the bytes at in[i].
*/
int lzp(unsigned char *in, int in_len, unsigned char *out) {
    int out_len = 0, i;
    int ht[1<<HASH_LEN] = {0};
    int hmask = (1<<HASH_LEN)-1;

    int h = 0;
    for (i = 0; i < in_len; i++) {
        if (ht[h] > 0) {
            int ml = match_len(in, i, in_len, &in[ht[h]]);
	    if (ml > 65535) ml = 65535;
            if (ml >= MIN_LEN) {
		if (ml <= 255) {
		    out[out_len++] = MATCH_CHAR;
		    out[out_len++] = ml;
		} else {
		    out[out_len++] = MATCH_CHAR+1;
		    out[out_len++] = ml>>8;
		    out[out_len++] = ml;
		}
		// Not as accurate, but a speed up on long matches
		// It's only around 5-10% total gain though so maybe
		// not worth the complexity (of spec).
#ifdef FAST_MODE
		if (ml > 4) {
		    i += ml-4;
		    ml = 4;
		}
#endif
                do {
                    ht[h] = i;
                    h = update_hash(h, in[i]);
                    i++;
                } while (--ml > 0);
                i--;
            } else {
                if (in[i] == MATCH_CHAR || in[i] == MATCH_CHAR+1) {
                    // encode as explicit zero length match instead
                    out[out_len++] = MATCH_CHAR;
                    out[out_len++] = 0;
                }
                out[out_len++] = in[i];
                ht[h] = i;
                h = update_hash(h, in[i]);
            }
        } else {
            out[out_len++] = in[i];
            ht[h] = i;
            h = update_hash(h, in[i]);
        }
    }

    return out_len;
}

int unlzp(unsigned char *in, int in_len, unsigned char *out) {
    int i, j;
    int ht[1<<HASH_LEN] = {0};
    int hmask = (1<<HASH_LEN)-1;

    int h = 0;
    for (i = j = 0; i < in_len; i++) {
	if (ht[h] > 0) {
	    int is_match = in[i++], ml = 0;
	    if (is_match == MATCH_CHAR) {
		ml = in[i++];
	    } else if (is_match == MATCH_CHAR+1) {
		ml = in[i++]<<8;
		ml += in[i++];
	    }
	    if (ml) {
		if (ht[h]+ml < j) {
		    memcpy(&out[j], &out[ht[h]], ml);
		} else {
		    int z;
		    for (z = 0; z < ml; z++)
			out[j+z] = out[ht[h]+z];
		}
#ifdef FAST_MODE
		if (ml > 4) {
		    j += ml-4;
		    ml = 4;
		}
#endif
 		do {
		    ht[h] = j;
		    h = update_hash(h, out[j]);
		    j++;
		} while (--ml > 0);
		i--;
	    } else {
		i -= (is_match != MATCH_CHAR && is_match != MATCH_CHAR+1);
		out[j] = in[i];
		ht[h] = j++;
		h = update_hash(h, in[i]);
	    }
	} else {
	    out[j] = in[i];
	    ht[h] = j++;
	    h = update_hash(h, in[i]);
	}
    }

    return j;
}
