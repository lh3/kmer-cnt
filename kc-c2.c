#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include "ketopt.h" // command-line argument parser

#include "kseq.h" // FASTA/Q parser
KSEQ_INIT(gzFile, gzread)

#include "khashl.h" // hash table
#define KC_BITS 10
#define kc_c2_eq(a, b) ((a)>>KC_BITS == (b)>>KC_BITS) // lower 8 bits for counts; higher bits for k-mer
#define kc_c2_hash(a) ((a)>>KC_BITS)
KHASHL_SET_INIT(, kc_c2_t, kc_c2, uint64_t, kc_c2_hash, kc_c2_eq)

#define CALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))

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
	int p;
	kc_c2_t **h;
} kc_c2x_t;

static kc_c2x_t *c2x_init(int p)
{
	int i;
	kc_c2x_t *h;
	CALLOC(h, 1);
	CALLOC(h->h, 1<<p);
	h->p = p;
	for (i = 0; i < 1<<p; ++i)
		h->h[i] = kc_c2_init();
	return h;
}

static inline void c2x_insert(kc_c2x_t *h, uint64_t y) // insert a k-mer $y to hash table $h
{
	int absent, pre = y & ((1<<h->p) - 1);
	kc_c2_t *g = h->h[pre];
	khint_t k;
	k = kc_c2_put(g, y>>h->p<<KC_BITS, &absent);
	if ((kh_key(g, k)&0xff) < 255) ++kh_key(g, k); // count if not saturated
}

static void count_seq(kc_c2x_t *h, int k, int len, char *seq) // insert k-mers in $seq to hash table $h
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
				c2x_insert(h, hash64(y, mask));
			}
		} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
	}
}

static kc_c2x_t *count_file(const char *fn, int k, int p)
{
	gzFile fp;
	kseq_t *ks;
	kc_c2x_t *h;
	if ((fp = gzopen(fn, "r")) == 0) return 0;
	ks = kseq_init(fp);
	h = c2x_init(p);
	while (kseq_read(ks) >= 0)
		count_seq(h, k, ks->seq.l, ks->seq.s);
	kseq_destroy(ks);
	gzclose(fp);
	return h;
}

static void print_hist(const kc_c2x_t *h)
{
	khint_t k;
	uint64_t cnt[256];
	int i;
	for (i = 0; i < 256; ++i) cnt[i] = 0;
	for (i = 0; i < 1<<h->p; ++i) {
		kc_c2_t *g = h->h[i];
		for (k = 0; k < kh_end(g); ++k)
			if (kh_exist(g, k))
				++cnt[kh_key(g, k)&0xff];
	}
	for (i = 1; i < 256; ++i)
		printf("%d\t%ld\n", i, (long)cnt[i]);
}

int main(int argc, char *argv[])
{
	kc_c2x_t *h;
	int i, c, k = 31, p = KC_BITS;
	ketopt_t o = KETOPT_INIT;
	while ((c = ketopt(&o, argc, argv, 1, "k:p:", 0)) >= 0)
		if (c == 'k') k = atoi(o.arg);
		else if (c == 'p') p = atoi(o.arg);
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: kc-c2 [-k %d] [-p %d] <in.fa>\n", k, p);
		return 1;
	}
	if (p < KC_BITS) {
		fprintf(stderr, "ERROR: -p should be at least %d\n", KC_BITS);
		return 1;
	}
	h = count_file(argv[o.ind], k, p);
	print_hist(h);
	for (i = 0; i < 1<<p; ++i)
		kc_c2_destroy(h->h[i]);
	free(h->h); free(h);
	return 0;
}
