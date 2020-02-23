#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include "ketopt.h" // command-line argument parser
#include "kthread.h" // multi-threading models: pipeline and multi-threaded for loop

#include "kseq.h" // FASTA/Q parser
KSEQ_INIT(gzFile, gzread)

#include "khashl.h" // hash table
#define kc_c4_eq(a, b) ((a)>>8 == (b)>>8) // lower 8 bits for counts; higher bits for k-mer
#define kc_c4_hash(a) ((a)>>8)
KHASHL_SET_INIT(, kc_c4_t, kc_c4, uint64_t, kc_c4_hash, kc_c4_eq)

#define CALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define MALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))
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
	kc_c4_t **h; // 1<<p hash tables
} kc_c4x_t;

static kc_c4x_t *c4x_init(int p)
{
	int i;
	kc_c4x_t *h;
	CALLOC(h, 1);
	MALLOC(h->h, 1<<p);
	h->p = p;
	for (i = 0; i < 1<<p; ++i)
		h->h[i] = kc_c4_init();
	return h;
}

typedef struct {
	int n, m;
	uint64_t *a;
} buf_c4_t;

static inline void c4x_insert_buf(buf_c4_t *buf, int p, uint64_t y) // insert a k-mer $y to a linear buffer
{
	int pre = y & ((1<<p) - 1);
	buf_c4_t *b = &buf[pre];
	if (b->n == b->m) {
		b->m = b->m < 8? 8 : b->m + (b->m>>1);
		REALLOC(b->a, b->m);
	}
	b->a[b->n++] = y;
}

static void count_seq_buf(buf_c4_t *buf, int k, int p, int len, const char *seq) // insert k-mers in $seq to linear buffer $buf
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
				c4x_insert_buf(buf, p, hash64(y, mask));
			}
		} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
	}
}

typedef struct { // global data structure for kt_pipeline()
	int k, block_len, n_thread;
	kseq_t *ks;
	kc_c4x_t *h;
} pldat_t;

typedef struct { // data structure for each step in kt_pipeline()
	pldat_t *p;
	int n, m, sum_len, nk;
	int *len;
	char **seq;
	buf_c4_t *buf;
} stepdat_t;

static void worker_for(void *data, long i, int tid) // callback for kt_for()
{
	stepdat_t *s = (stepdat_t*)data;
	buf_c4_t *b = &s->buf[i];
	kc_c4_t *h = s->p->h->h[i];
	int j, p = s->p->h->p;
	for (j = 0; j < b->n; ++j) {
		khint_t k;
		int absent;
		k = kc_c4_put(h, b->a[j]>>p<<8, &absent);
		if ((kh_key(h, k)&0xff) < 255) ++kh_key(h, k);
	}
}

static void *worker_pipeline(void *data, int step, void *in) // callback for kt_pipeline()
{
	pldat_t *p = (pldat_t*)data;
	if (step == 0) { // step 1: read a block of sequences
		int ret;
		stepdat_t *s;
		CALLOC(s, 1);
		s->p = p;
		while ((ret = kseq_read(p->ks)) >= 0) {
			int l = p->ks->seq.l;
			if (l < p->k) continue;
			if (s->n == s->m) {
				s->m = s->m < 16? 16 : s->m + (s->n>>1);
				REALLOC(s->len, s->m);
				REALLOC(s->seq, s->m);
			}
			MALLOC(s->seq[s->n], l);
			memcpy(s->seq[s->n], p->ks->seq.s, l);
			s->len[s->n++] = l;
			s->sum_len += l;
			s->nk += l - p->k + 1;
			if (s->sum_len >= p->block_len)
				break;
		}
		if (s->sum_len == 0) free(s);
		else return s;
	} else if (step == 1) { // step 2: extract k-mers
		stepdat_t *s = (stepdat_t*)in;
		int i, n = 1<<p->h->p, m;
		CALLOC(s->buf, n);
		m = (int)(s->nk * 1.2 / n) + 1;
		for (i = 0; i < n; ++i) {
			s->buf[i].m = m;
			MALLOC(s->buf[i].a, m);
		}
		for (i = 0; i < s->n; ++i) {
			count_seq_buf(s->buf, p->k, p->h->p, s->len[i], s->seq[i]);
			free(s->seq[i]);
		}
		free(s->seq); free(s->len);
		return s;
	} else if (step == 2) { // step 3: insert k-mers to hash table
		stepdat_t *s = (stepdat_t*)in;
		int i, n = 1<<p->h->p;
		kt_for(p->n_thread, worker_for, s, n);
		for (i = 0; i < n; ++i) free(s->buf[i].a);
		free(s->buf); free(s);
	}
	return 0;
}

static kc_c4x_t *count_file(const char *fn, int k, int p, int block_size, int n_thread)
{
	pldat_t pl;
	gzFile fp;
	if ((fp = gzopen(fn, "r")) == 0) return 0;
	pl.ks = kseq_init(fp);
	pl.k = k;
	pl.n_thread = n_thread;
	pl.h = c4x_init(p);
	pl.block_len = block_size;
	kt_pipeline(3, worker_pipeline, &pl, 3);
	kseq_destroy(pl.ks);
	gzclose(fp);
	return pl.h;
}

typedef struct {
	uint64_t c[256];
} buf_cnt_t;

typedef struct {
	const kc_c4x_t *h;
	buf_cnt_t *cnt;
} hist_aux_t;

static void worker_hist(void *data, long i, int tid) // callback for kt_for()
{
	hist_aux_t *a = (hist_aux_t*)data;
	uint64_t *cnt = a->cnt[tid].c;
	kc_c4_t *g = a->h->h[i];
	khint_t k;
	for (k = 0; k < kh_end(g); ++k)
		if (kh_exist(g, k))
			++cnt[kh_key(g, k)&0xff];
}

static void print_hist(const kc_c4x_t *h, int n_thread)
{
	hist_aux_t a;
	uint64_t cnt[256];
	int i, j;
	a.h = h;
	CALLOC(a.cnt, n_thread);
	kt_for(n_thread, worker_hist, &a, 1<<h->p);
	for (i = 0; i < 256; ++i) cnt[i] = 0;
	for (j = 0; j < n_thread; ++j)
		for (i = 0; i < 256; ++i)
			cnt[i] += a.cnt[j].c[i];
	free(a.cnt);
	for (i = 1; i < 256; ++i)
		printf("%d\t%ld\n", i, (long)cnt[i]);
}

int main(int argc, char *argv[])
{
	kc_c4x_t *h;
	int i, c, k = 31, p = 8, block_size = 10000000, n_thread = 4;
	ketopt_t o = KETOPT_INIT;
	while ((c = ketopt(&o, argc, argv, 1, "k:p:b:t:", 0)) >= 0) {
		if (c == 'k') k = atoi(o.arg);
		else if (c == 'p') p = atoi(o.arg);
		else if (c == 'b') block_size = atoi(o.arg);
		else if (c == 't') n_thread = atoi(o.arg);
	}
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: kc-c4 [options] <in.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT     k-mer size [%d]\n", k);
		fprintf(stderr, "  -p INT     prefix length [%d]\n", p);
		fprintf(stderr, "  -b INT     block size [%d]\n", block_size);
		fprintf(stderr, "  -t INT     number of worker threads [%d]\n", n_thread);
		return 1;
	}
	if (p < 8) {
		fprintf(stderr, "ERROR: -p should be at least 8\n");
		return 1;
	}
	h = count_file(argv[o.ind], k, p, block_size, n_thread);
	print_hist(h, n_thread);
	for (i = 0; i < 1<<p; ++i) kc_c4_destroy(h->h[i]);
	free(h->h); free(h);
	return 0;
}
