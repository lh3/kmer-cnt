#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include "ketopt.h" // command-line argument parser

#include "kseq.h" // FASTA/Q parser
KSEQ_INIT(gzFile, gzread)

#include "khashl.h" // hash table
KHASHL_MAP_INIT(, kc_c1_t, kc_c1, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)

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

static void count_seq(kc_c1_t *h, int k, int len, char *seq) // insert k-mers in $seq to hash table $h
{
	int i, l;
	uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
		int absent, c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			if (++l >= k) { // we find a k-mer
				khint_t itr;
				uint64_t y = x[0] < x[1]? x[0] : x[1];
				itr = kc_c1_put(h, y, &absent); // only add one strand!
				if (absent) kh_val(h, itr) = 0;
				++kh_val(h, itr);
			}
		} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
	}
}

static kc_c1_t *count_file(const char *fn, int k)
{
	gzFile fp;
	kseq_t *ks;
	kc_c1_t *h;
	if ((fp = gzopen(fn, "r")) == 0) return 0;
	ks = kseq_init(fp);
	h = kc_c1_init();
	while (kseq_read(ks) >= 0)
		count_seq(h, k, ks->seq.l, ks->seq.s);
	kseq_destroy(ks);
	gzclose(fp);
	return h;
}

static void print_hist(const kc_c1_t *h)
{
	khint_t k;
	uint64_t cnt[256];
	int i;
	for (i = 0; i < 256; ++i) cnt[i] = 0;
	for (k = 0; k < kh_end(h); ++k)
		if (kh_exist(h, k))
			++cnt[kh_val(h, k) < 256? kh_val(h, k) : 255];
	for (i = 1; i < 256; ++i)
		printf("%d\t%ld\n", i, (long)cnt[i]);
}

int main(int argc, char *argv[])
{
	kc_c1_t *h;
	int c, k = 31;
	ketopt_t o = KETOPT_INIT;
	while ((c = ketopt(&o, argc, argv, 1, "k:", 0)) >= 0)
		if (c == 'k') k = atoi(o.arg);
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: kc-c1 [-k %d] <in.fa>\n", k);
		return 1;
	}
	h = count_file(argv[o.ind], k);
	print_hist(h);
	kc_c1_destroy(h);
	return 0;
}
