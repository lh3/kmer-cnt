#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include "ketopt.h" // command-line argument parser
#include "kthread.h"

#include "kseq.h" // FASTA/Q parser
KSEQ_INIT(gzFile, gzread)

#include "khashl.h" // hash table
#define KC_BITS 10
#define kc_c3_eq(a, b) ((a)>>KC_BITS == (b)>>KC_BITS) // lower 8 bits for counts; higher bits for k-mer
#define kc_c3_hash(a) ((a)>>KC_BITS)
KHASHL_SET_INIT(, kc_c3_t, kc_c3, uint64_t, kc_c3_hash, kc_c3_eq)

#define CALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))

const unsigned char seq_nt4_table[256] = { // translate ACGT to 0123
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static inline uint64_t hash64(uint64_t key, uint64_t mask) // invertible integer hash function
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

typedef struct {
	int p; // suffix length; at least 8
	kc_c3_t **h; // 1<<p hash tables
} kc_c3x_t;

static kc_c3x_t *c3x_init(int p)
{
	int i;
	kc_c3x_t *h;
	CALLOC(h, 1);
	CALLOC(h->h, 1<<p);
	h->p = p;
	for (i = 0; i < 1<<p; ++i)
		h->h[i] = kc_c3_init();
	return h;
}

static inline void c3x_insert(kc_c3x_t *h, uint64_t y) // insert a k-mer $y to hash table $h
{
	int absent, pre = y & ((1<<h->p) - 1);
	kc_c3_t *g = h->h[pre];
	khint_t k;
	k = kc_c3_put(g, y>>h->p<<KC_BITS, &absent);
	if ((kh_key(g, k)&0xff) < 255) ++kh_key(g, k); // count if not saturated
}

static void count_seq(kc_c3x_t *h, int k, int len, char *seq) // insert k-mers in $seq to hash table $h
{
	int i, l;
	uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			if (++l >= k) { // we find a k-mer
				uint64_t y = x[0] < x[1]? x[0] : x[1];
				c3x_insert(h, hash64(y, mask));
			}
		} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
	}
}

typedef struct {
	int k, block_len;
	kseq_t *ks;
	kc_c3x_t *h;
} pldat_t;

typedef struct {
	int n, m, sum_len;
	int *len;
	char **seq;
} stepdat_t;

static void *worker_pipeline(void *data, int step, void *in) // callback for kt_pipeline()
{
	pldat_t *p = (pldat_t*)data;
	if (step == 0) { // step 1: read a block of sequences
		int ret;
		stepdat_t *s;
		CALLOC(s, 1);
		while ((ret = kseq_read(p->ks)) >= 0) {
			int l = p->ks->seq.l;
			if (s->n == s->m) {
				s->m = s->m < 16? 16 : s->m + (s->n>>1);
				REALLOC(s->len, s->m);
				REALLOC(s->seq, s->m);
			}
			CALLOC(s->seq[s->n], l);
			memcpy(s->seq[s->n], p->ks->seq.s, l);
			s->len[s->n++] = l;
			s->sum_len += l;
			if (s->sum_len >= p->block_len)
				break;
		}
		if (s->sum_len == 0) free(s);
		else return s;
	} else if (step == 1) { // step 2: count k-mers
		stepdat_t *s = (stepdat_t*)in;
		int i;
		for (i = 0; i < s->n; ++i) {
			count_seq(p->h, p->k, s->len[i], s->seq[i]);
			free(s->seq[i]);
		}
		free(s->seq); free(s->len); free(s);
	}
	return 0;
}

static kc_c3x_t *count_file(const char *fn, int k, int p, int block_size)
{
	pldat_t pl;
	gzFile fp;
	if ((fp = gzopen(fn, "r")) == 0) return 0;
	pl.ks = kseq_init(fp);
	pl.k = k;
	pl.h = c3x_init(p);
	pl.block_len = block_size;
	kt_pipeline(2, worker_pipeline, &pl, 2);
	kseq_destroy(pl.ks);
	gzclose(fp);
	return pl.h;
}

static void print_hist(const kc_c3x_t *h)
{
	khint_t k;
	uint64_t cnt[256];
	int i;
	for (i = 0; i < 256; ++i) cnt[i] = 0;
	for (i = 0; i < 1<<h->p; ++i) {
		kc_c3_t *g = h->h[i];
		for (k = 0; k < kh_end(g); ++k)
			if (kh_exist(g, k))
				++cnt[kh_key(g, k)&0xff];
	}
	for (i = 1; i < 256; ++i)
		printf("%d\t%ld\n", i, (long)cnt[i]);
}

int main(int argc, char *argv[])
{
	kc_c3x_t *h;
	int i, c, k = 31, p = KC_BITS, block_size = 10000000;
	ketopt_t o = KETOPT_INIT;
	while ((c = ketopt(&o, argc, argv, 1, "k:p:b:", 0)) >= 0) {
		if (c == 'k') k = atoi(o.arg);
		else if (c == 'p') p = atoi(o.arg);
		else if (c == 'b') block_size = atoi(o.arg);
	}
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: kc-c3 [-k %d] [-p %d] [-b %d] <in.fa>\n", k, p, block_size);
		return 1;
	}
	if (p < KC_BITS) {
		fprintf(stderr, "ERROR: -p should be at least %d\n", KC_BITS);
		return 1;
	}
	h = count_file(argv[o.ind], k, p, block_size);
	print_hist(h);
	for (i = 0; i < 1<<p; ++i)
		kc_c3_destroy(h->h[i]);
	free(h->h); free(h);
	return 0;
}
