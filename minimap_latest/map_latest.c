#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <pthread.h>
#include "bseq.h"
#include "minimap.h"
#include "kvec.h"
#include "sdust.h"

FILE *fpUnMapped;

const unsigned char seq_nt16_table[256] = {
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
     1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0 /*=*/,15,15,
    15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
    15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,

    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};
void mm_mapopt_init(mm_mapopt_t *opt)
{
	opt->radius = 100;//vasu 500
	opt->max_gap = 10000;
	opt->min_cnt = 4;
	opt->min_match = 40;
	opt->sdust_thres = 0;
	opt->flag = MM_F_WITH_REP;
	opt->merge_frac = 1;//vasu .5
}

/****************************
 * Find approxiate mappings *
 ****************************/

struct mm_tbuf_s { // per-thread buffer
	mm128_v mini; // query minimizers
	mm128_v coef; // Hough transform coefficient
	mm128_v intv; // intervals on sorted coef
	uint32_v reg2mini;
	uint32_v rep_aux;
	sdust_buf_t *sdb;
	// the following are for computing LIS
	uint32_t n, m;
	uint64_t *a;
	size_t *b, *p;
	// final output
	kvec_t(mm_reg1_t) reg;
};

mm_tbuf_t *mm_tbuf_init()
{
	mm_tbuf_t *b;
	b = (mm_tbuf_t*)calloc(1, sizeof(mm_tbuf_t));
	b->sdb = sdust_buf_init();
	return b;
}

void mm_tbuf_destroy(mm_tbuf_t *b)
{
	if (b == 0) return;
	free(b->mini.a); free(b->coef.a); free(b->intv.a); free(b->reg.a); free(b->reg2mini.a); free(b->rep_aux.a);
	free(b->a); free(b->b); free(b->p);
	sdust_buf_destroy(b->sdb);
	free(b);
}

#include "ksort.h"
#define sort_key_64(a) (a)
KRADIX_SORT_INIT(64, uint64_t, sort_key_64, 8) 
#define lt_low32(a, b) ((uint32_t)(a) < (uint32_t)(b))
KSORT_INIT(low32lt, uint64_t, lt_low32)
#define gt_low32(a, b) ((uint32_t)(a) > (uint32_t)(b))
KSORT_INIT(low32gt, uint64_t, gt_low32)

/* TODO: drop_rep() is not robust. For all-vs-all mapping but without the -S
 * flag, all minimizers have at least one hit. The _thres_ computed below will
 * be highly skewed. Some improvements need to be made. */

static void drop_rep(mm_tbuf_t *b, int min_cnt)
{
	int i, j, n, m;
	uint32_t thres;
	b->rep_aux.n = 0;
	for (i = 0; i < b->mini.n; ++i)
		if (b->mini.a[i].y>>32)
			kv_push(uint32_t, b->rep_aux, b->mini.a[i].y>>32);
	if (b->rep_aux.n < 3) return;
	thres = (uint32_t)(ks_ksmall_uint32_t(b->rep_aux.n, b->rep_aux.a, b->rep_aux.n>>1) * MM_DEREP_Q50 + .499);
	for (i = n = m = 0; i < b->reg.n; ++i) {
		int cnt = 0, all_cnt = b->reg.a[i].cnt;
		for (j = 0; j < all_cnt; ++j)
			if (b->mini.a[b->reg2mini.a[m + j]].y>>32 <= thres)
				++cnt;
		if (cnt >= min_cnt)
			b->reg.a[n++] = b->reg.a[i];
		m += all_cnt;
	}
//	printf("%ld=>%d\t%d\n", b->reg.n, n, thres);
	b->reg.n = n;
}

static void proc_intv(mm_tbuf_t *b, int which, int k, int min_cnt, int max_gap)
{
	int i, j, l_lis, rid = -1, rev = 0, start = b->intv.a[which].y, end = start + b->intv.a[which].x;

	// make room for arrays needed by LIS (longest increasing sequence)
	if (end - start > b->m) {
		b->m = end - start;
		kv_roundup32(b->m);
		b->a = (uint64_t*)realloc(b->a, b->m * 8);
		b->b = (size_t*)realloc(b->b, b->m * sizeof(size_t));
		b->p = (size_t*)realloc(b->p, b->m * sizeof(size_t));
	}

	// prepare the input array _a_ for LIS
	b->n = 0;
	for (i = start; i < end; ++i)
		if (b->coef.a[i].x != UINT64_MAX)
			b->a[b->n++] = b->coef.a[i].y, rid = b->coef.a[i].x << 1 >> 33, rev = b->coef.a[i].x >> 63;
	if (b->n < min_cnt) return;
	radix_sort_64(b->a, b->a + b->n);

	// find the longest increasing sequence
	l_lis = rev? ks_lis_low32gt(b->n, b->a, b->b, b->p) : ks_lis_low32lt(b->n, b->a, b->b, b->p); // LIS
	if (l_lis < min_cnt) return;
	for (i = 1, j = 1; i < l_lis; ++i) // squeeze out minimizaers reused in the LIS sequence
		if (b->a[b->b[i]]>>32 != b->a[b->b[i-1]]>>32)
			b->a[b->b[j++]] = b->a[b->b[i]];
	l_lis = j;
	if (l_lis < min_cnt) return;

	// convert LISes to regions; possibly break an LIS at a long gaps
	for (i = 1, start = 0; i <= l_lis; ++i) {
		int32_t qgap = i == l_lis? 0 : ((uint32_t)b->mini.a[b->a[b->b[i]]>>32].y>>1) - ((uint32_t)b->mini.a[b->a[b->b[i-1]]>>32].y>>1);
		if (i == l_lis || (qgap > max_gap && abs((int32_t)b->a[b->b[i]] - (int32_t)b->a[b->b[i-1]]) > max_gap)) {
			if (i - start >= min_cnt) {
				uint32_t lq = 0, lr = 0, eq = 0, er = 0, sq = 0, sr = 0;
				mm_reg1_t *r;
				kv_pushp(mm_reg1_t, b->reg, &r);
				r->rid = rid, r->rev = rev, r->cnt = i - start, r->rep = 0;
				r->qs = ((uint32_t)b->mini.a[b->a[b->b[start]]>>32].y>>1) - (k - 1);
				r->qe = ((uint32_t)b->mini.a[b->a[b->b[i-1]]>>32].y>>1) + 1;
				r->rs = rev? (uint32_t)b->a[b->b[i-1]] : (uint32_t)b->a[b->b[start]];
				r->re = rev? (uint32_t)b->a[b->b[start]] : (uint32_t)b->a[b->b[i-1]];
				r->rs -= k - 1;
				r->re += 1;
				for (j = start; j < i; ++j) { // count the number of times each minimizer is used
					int jj = b->a[b->b[j]]>>32;
					b->mini.a[jj].y += 1ULL<<32;
					kv_push(uint32_t, b->reg2mini, jj); // keep minimizer<=>reg mapping for derep
				}
				for (j = start; j < i; ++j) { // compute ->len
					uint32_t q = ((uint32_t)b->mini.a[b->a[b->b[j]]>>32].y>>1) - (k - 1);
					uint32_t r = (uint32_t)b->a[b->b[j]];
					r = !rev? r - (k - 1) : (0x80000000U - r);
					if (r > er) lr += er - sr, sr = r, er = sr + k;
					else er = r + k;
					if (q > eq) lq += eq - sq, sq = q, eq = sq + k;
					else eq = q + k;
				}
				lr += er - sr, lq += eq - sq;
				r->len = lr < lq? lr : lq;
			}
			start = i;
		}
	}
}

// merge or add a Hough interval; only used by get_reg()
static inline void push_intv(mm128_v *intv, int start, int end, float merge_frac)
{
	mm128_t *p;
	if (intv->n > 0) { // test overlap
		int last_start, last_end, min;
		p = &intv->a[intv->n-1];
		last_start = p->y, last_end = p->x + last_start;
		min = end - start < last_end - last_start? end - start : last_end - last_start;
		if (last_end > start && last_end - start > min * merge_frac) { // large overlap; then merge //vasu: dont merge intervals
			p->x = end - last_start;
//			printf ("Dont push\n");
			return;
		}
	}
	kv_pushp(mm128_t, *intv, &p); // a new interval
	p->x = end - start, p->y = start;
}

// find mapping regions from a list of minimizer hits
static void get_reg(mm_tbuf_t *b, int radius, int k, int min_cnt, int max_gap, float merge_frac, int flag)
{
	const uint64_t v_kept = ~(1ULL<<31), v_dropped = 1ULL<<31;
	mm128_v *c = &b->coef;
	int i, j, start = 0, iso_dist = radius * 2;
//	printf ("c->n %lu\n", c->n);
	if (c->n < min_cnt) return;

	// drop isolated minimizer hits
	if (flag&MM_F_NO_ISO) {
//		printf ("get_reg: MM_F_NO_ISO\n");
		for (i = 0; i < c->n; ++i) c->a[i].y |= v_dropped;
		for (i = 1; i < c->n; ++i) {
			uint64_t x = c->a[i].x;
			int32_t rpos = (uint32_t)c->a[i].y;
			for (j = i - 1; j >= 0 && x - c->a[j].x < radius; --j) {
				int32_t y = c->a[j].y;
				if (abs(y - rpos) < iso_dist) {
					c->a[i].y &= v_kept, c->a[j].y &= v_kept;
					break;
				}
			}
		}
		for (i = j = 0; i < c->n; ++i) // squeeze out hits still marked as v_dropped
			if ((c->a[i].y&v_dropped) == 0)
				c->a[j++] = c->a[i];
		c->n = j;
	}
//	printf ("get_reg: out of if\n");
	// identify (possibly overlapping) intervals within _radius_; an interval is a cluster of hits
	b->intv.n = 0;
	for (i = 1; i < c->n; ++i) {
//		printf ("i %d x %llx y %llx, start %d, posDiff %lu radius %d\n", i, c->a[i].x, c->a[i].y, start, c->a[i].x - c->a[start].x, radius);
		if (c->a[i].x - c->a[start].x > radius) {
			if (i - start >= min_cnt) 
			{
				push_intv(&b->intv, start, i, merge_frac);
//				printf ("\tb->intv.n %lu numMinimizerHits %d\n", b->intv.n, i-start);
			}
			for (++start; start < i && c->a[i].x - c->a[start].x > radius; ++start);
		}
	}
	if (i - start >= min_cnt) push_intv(&b->intv, start, i, merge_frac);

	// sort by the size of the interval
	radix_sort_128x(b->intv.a, b->intv.a + b->intv.n);

	// generate hits, starting from the largest interval
	b->reg2mini.n = 0;
	int numCandidates=0;
	for (i = b->intv.n - 1; i >= 0; --i) 
	{
//		if (numCandidates == 10) break;
		proc_intv(b, i, k, min_cnt, max_gap);
		numCandidates++;
	}

	// post repeat removal
	if (!(flag&MM_F_WITH_REP)) drop_rep(b, min_cnt);
}

const mm_reg1_t *mm_map(const mm_idx_t *mi, int l_seq, const char *seq, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *name)
{
	int j, n_dreg = 0, u = 0;
	const uint64_t *dreg = 0;

	b->mini.n = b->coef.n = 0;
	mm_sketch(seq, l_seq, mi->w, mi->k, 0, &b->mini);
//	printf ("mm_map name %s, w %d, k %d, numMinimizers %lu, opt->sdust_thres %d\n", name, mi->w, mi->k, b->mini.n, opt->sdust_thres);
	if (opt->sdust_thres > 0)
		dreg = sdust_core((const uint8_t*)seq, l_seq, opt->sdust_thres, 64, &n_dreg, b->sdb);
	for (j = 0; j < b->mini.n; ++j) {
		int k, n;
		const uint64_t *r;
		int32_t qpos = (uint32_t)b->mini.a[j].y>>1, strand = b->mini.a[j].y&1;
//		printf ("j %d, qpos %d\n", j, qpos);
		b->mini.a[j].y = b->mini.a[j].y<<32>>32; // clear the rid field
		if (dreg && n_dreg) { // test complexity
			int s = qpos - (mi->k - 1), e = s + mi->k;
			while (u < n_dreg && (uint32_t)dreg[u] <= s) ++u;
			if (u < n_dreg && dreg[u]>>32 < e) {
				int v, l = 0;
				for (v = u; v < n_dreg && dreg[v]>>32 < e; ++v) { // iterate over LCRs overlapping this minimizer
					int ss = s > dreg[v]>>32? s : dreg[v]>>32;
					int ee = e < (uint32_t)dreg[v]? e : (uint32_t)dreg[v];
					l += ee - ss;
				}
				if (l > mi->k>>1) continue;
			}
		}
		r = mm_idx_get(mi, b->mini.a[j].x, &n);
//		printf ("MINI %lx,%lu n %d max_occ %d\n", b->mini.a[j].x, b->mini.a[j].y>>1, n, mi->max_occ);
		if (n > mi->max_occ) continue;//vasu
		for (k = 0; k < n; ++k) {
			int32_t rpos = (uint32_t)r[k] >> 1;
			mm128_t *p;
			if (name && (opt->flag&MM_F_NO_SELF) && mi->name && strcmp(name, mi->name[r[k]>>32]) == 0 && rpos == qpos)
				continue;
			if (name && (opt->flag&MM_F_AVA) && mi->name && strcmp(name, mi->name[r[k]>>32]) > 0)
				continue;
			kv_pushp(mm128_t, b->coef, &p);
			if ((r[k]&1) == strand) { // forward strand
				p->x = (uint64_t)r[k] >> 32 << 32 | (0x80000000U + rpos - qpos);
				p->y = (uint64_t)j << 32 | rpos;
//				printf ("1-hash %llx - J %d k %d X %llx (%lx-%lx) Y %llx(%lx-%lx)\n", r[k], j, k, p->x, (uint64_t)r[k] >> 32, (0x80000000U + rpos - qpos), p->y, j, rpos);
			} else { // reverse strand
				p->x = (uint64_t)r[k] >> 32 << 32 | (rpos + qpos) | 1ULL<<63;
				p->y = (uint64_t)j << 32 | rpos;
//				printf ("2-hash %llx - J %d k %d X %llx (%lx-%lx) Y %llx(%lx-%lx)\n", r[k], j, k, p->x, (uint64_t)r[k] >> 32, (rpos + qpos), p->y, j, rpos);
			}
		}
	}
	radix_sort_128x(b->coef.a, b->coef.a + b->coef.n);
	b->reg.n = 0;
//	printf ("b->coef.n %lu opt->min_cnt %d, opt->radius %d\n", b->coef.n, opt->min_cnt, opt->radius);
	get_reg(b, opt->radius, mi->k, opt->min_cnt, opt->max_gap, opt->merge_frac, opt->flag);
//	printf ("b->reg.n %lu\n", b->reg.n);
	*n_regs = b->reg.n;
	return b->reg.a;
}

/**************************
 * Multi-threaded mapping *
 **************************/

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

typedef struct pipeline_t{
	int batch_size, n_processed, n_threads;
	const mm_mapopt_t *opt;
	bseq_file_t **fp1;
	bseq_file_t **fp2;
	const mm_idx_t *mi;
	char **refFasta;
	char *oprefix;
	int numOutChunks;
	int maxChunkSize;
	int numMaxChunks;
	FILE **aebFp;
	int *sizeToRead, *sizeToRead2;
	vec_fullRec *vfr;
	char **unmapped;
	size_t *u_m, *u_n;
	char **fastqBuffer1;
	size_t *fastqSize1;
	char **fastqBuffer2;
	size_t *fastqSize2;
} pipeline_t;

typedef struct {
	pipeline_t *p;
    int *n_seq;
	bseq1_t *seq1;
	bseq1_t *seq2;
	bseq1_t **seqs1;
	bseq1_t **seqs2;

	int *n_reg;
	int **n_regs;
	mm_reg1_t **reg;
	mm_reg1_t ***regs;

	mm_tbuf_t **buf1;
	mm_tbuf_t **buf2;
} step_t;

static int reg_compare (const void *p1, const void *q1)
{
	mm_reg1_t * p = (mm_reg1_t *)p1;
	mm_reg1_t * q = (mm_reg1_t *)q1;
	if ((p->score < q->score) || ((p->score == q->score) && (p->rid < q->rid)) || ((p->score == q->score) && (p->rid == q->rid) && (p->rs < q->rs)))
		return 1;
	else if ((p->score > q->score) || ((p->score == q->score) && (p->rid > q->rid)) || ((p->score == q->score) && (p->rid == q->rid) && (p->rs >  q->rs)))
		return -1;
	else
		return 0;
}
void assign_score (char **refFasta, bseq1_t *seq, mm_reg1_t *reg, int n_regs)
{
	int i1=0, j1=0;
	for (i1=0; i1<n_regs; i1++)
	{
		mm_reg1_t *pr = &(reg[i1]);
		pr->score=0;
		if (pr->rev == 0)//'+' strand
		{
			int rpos = pr->rs - pr->qs;
			for (j1=0; j1<seq->l_seq; j1++)
			{
				if (seq->seq[j1] > 96) seq->seq[j1] -= 32;
				if (rpos < 0) 
				{
					rpos++;
					continue;
				}
				if (seq->seq[j1] == refFasta[pr->rid][rpos]) pr->score++;
				rpos++;
			}
		}
		else
		{
			int rpos = pr->rs - (seq->l_seq - pr->qe);
			for (j1=(seq->l_seq-1); j1>=0; j1--)
			{
				if (seq->seq[j1] > 96) seq->seq[j1] -= 32;
				char rc = seq->seq[j1];
				if (rc != 'N')
					rc = "TGAC"[(rc >> 1) & 3];
				if (rpos < 0) 
				{
					rpos++;
					continue;
				}
				if (rc==refFasta[pr->rid][rpos]) pr->score++;
				rpos++;
			}
		}
		pr->score = seq->l_seq - 4*(seq->l_seq - pr->score);
		if (pr->score < 0) pr->score = 0;
	}
}

int calcMapQ (mm_reg1_t *pr, int pr_n, int sub_n, int i, int rLen)
{
	int mapQ=0;
	double identity = 1. - (double)(rLen - pr[0].score) / 5 / rLen;
	double tmp = (double)3 / 4.61512;
	tmp *= identity * identity;
//	printf ("identity %lf tmp %lf\n", identity, tmp);
//	int i=0;
//	while ((i< pr_n) && (pr[i].score == pr[0].score) && (pr[i].rs == pr[0].rs)) i++;
//	printf ("i1 %d/%d ", i, pr_n);
	if (i == pr_n)
		mapQ = 60;
	else if (pr[0].score == pr[i].score)
		mapQ = 0;
	else
	{
		int sub = (pr[i].score < 19) ? 19: pr[i].score;
//		printf ("Score %d sub %d sub_n %d\n", pr[0].score, sub, sub_n);
		mapQ = (int)(6.02 * (pr[0].score - pr[i].score) * tmp * tmp + .499);
//		printf ("tmp mapQ %d\n", mapQ);
		if (sub_n > 0)
			mapQ -= (int) (4.343 * log(sub_n  + 1) + .499);
//		printf ("final mapQ %d\n", mapQ);
		if (mapQ > 60) mapQ = 60;
		if (mapQ < 0) mapQ =0;
	}
	return mapQ;
}

int selectAlnRecs (mm_reg1_t *regs1, mm_reg1_t *regs2, int n_regs1, int n_regs2, int rLen1, int rLen2, int *mq1, int *mq2)
{
	int posDiff=100000;
	int i=1;
	while ((i<n_regs1) && (regs1[i].rs == regs1[0].rs)) i++;
	int sub1_n=0, sub2_n=0, j=0, qual1=0, qual2=0;
	for (j=0; j<n_regs1; j++)
	{
		if ((j > i) && (regs1[j-1].rs != regs1[j].rs) && ((regs1[0].score - regs1[j].score) <= 7)) sub1_n++;
	}
	qual1 = calcMapQ(regs1, n_regs1, sub1_n, i, rLen1);
	regs1->mapQ = qual1;
	*mq1 = qual1;

	if (rLen2 == 0) return -1;
	int i1=0;
	for (i=0; i<n_regs2; i++)
	{
		if ((regs2[i].rid == regs1[0].rid) && (abs (regs2[i].rs - regs1[0].rs) < posDiff) && (regs2[i].rev != regs1[0].rev) && (regs2[i].score >= (rLen2 -4*4)) && (posDiff > 300) && ((posDiff == 100000) || ((regs2[i1].score-regs2[i].score) <= 2*4)))
		{
			i1=i;
			posDiff = abs (regs2[i].rs - regs1[0].rs);
		}
	}
	if (posDiff < 100000)
	{
		if (i1==0)
		{
			int sub2Ind = 0, sub2PosDiff=100000;
			for (i=i1+1; i<n_regs2; i++)
			{
				if ((regs2[i-1].rs != regs2[i].rs) && (regs2[i].rid == regs1[0].rid) && (regs2[i].rev != regs1[0].rev) && (abs (abs (regs2[i].rs - regs1[0].rs) - abs(regs2[i1].rs - regs1[0].rs)) < sub2PosDiff) && ((sub2PosDiff == 100000) || (regs2[i].score >= (rLen2 -10*4))))
				{
					sub2Ind = i;
					sub2PosDiff = abs (abs (regs2[i].rs - regs1[0].rs) - abs(regs2[i1].rs - regs1[0].rs));
					if (regs2[i].score == regs2[i1].score) break;
				} 
			}
			if (sub2Ind == 0) sub2Ind = n_regs2;
			for (j=0; j<n_regs2; j++)
			{
				if ((j >= sub2Ind) && (regs2[j-1].rs != regs2[j].rs) && ((regs2[0].score - regs2[j].score) <= 7)) sub2_n++; 
			}
			qual2 = calcMapQ(regs2, n_regs2, sub2_n, sub2Ind, rLen2);
			regs2->mapQ = qual2;
			*mq2 = qual2;
		}
		else
		{
			regs2[i1].mapQ = 0;
			*mq2 = 0;
		}
		return i1;
	}
	else
		return -1;
}

void create_aeb (bseq1_t *seq, mm_reg1_t *r, fullRec *fR)
{
	fR->pos = r->rs - r->qs+1;
	fR->flag = 0;
	if (r->rev) 
	{
		fR->pos = r->rs - (seq->l_seq - r->qe) +1;
		fR->flag = 16;
	}
	fR->qual = r->mapQ;
	fR->matchLen = seq->l_seq;
	int i=0;
	if (!r->rev)
	{
		for (i=0; i<seq->l_seq; i++)
		{
			int seqInd = i/2;
			int ofs = i&1;
			if (ofs==0)
				fR->seq[seqInd]=(uint8_t)seq_nt16_table[seq->seq[i]]<<4;
			else
				fR->seq[seqInd] |= (uint8_t)seq_nt16_table[seq->seq[i]];
//			printf ("Contig %d, i %d , %c , seqInd %d - %x\n", r->rid, i, seq->seq[i], seqInd, fR->seq[seqInd]);
			fR->quals[i] = seq->qual[i] - 33;
		}
	}
	else
	{
		for (i=seq->l_seq-1; i>=0; i--)
		{
			char rc = seq->seq[i];
			if (rc != 'N')
				rc = "TGAC"[(rc >> 1) & 3];
			int seqInd = (seq->l_seq-1 - i)/2;
			int ofs = (seq->l_seq-1 - i) & 1;
			if (ofs==0)
				fR->seq[seqInd]=(uint8_t)seq_nt16_table[rc]<<4;
			else
				fR->seq[seqInd] |= (uint8_t)seq_nt16_table[rc];
			fR->quals[seq->l_seq-1 - i] = seq->qual[i] - 33;
		}
		
	}
}

static void worker_for(void *_data, long i, int tid) // kt_for() callback
{
//	int i1=0, j1=0;
    step_t *step = (step_t*)_data;
	char **refFasta = step->p->refFasta;
	int maxChunkSize = step->p->maxChunkSize;
	int numMaxChunks = step->p->numMaxChunks;
	int n_threads = step->p->n_threads;
	int numContigs = step->p->mi->n;
	vec_fullRec *vfr = step->p->vfr;
	vfr += tid*numContigs*numMaxChunks;

	mm_reg1_t *regs[2]={NULL,NULL};
	int n_regs1, n_regs2=0;

	regs[0] = mm_map(step->p->mi, step->seq1[i].l_seq, step->seq1[i].seq, &n_regs1, step->buf1[tid], step->p->opt, step->seq1[i].name);
	step->n_reg[2*i] = n_regs1;
	if (step->seq2 != NULL)
	{
		regs[1] = mm_map(step->p->mi, step->seq2[i].l_seq, step->seq2[i].seq, &n_regs2, step->buf2[tid], step->p->opt, step->seq2[i].name);
		step->n_reg[2*i+1] = n_regs2;
	}
	
//	printf ("tid %d worker_for name %s, n_regs1 %d, n_regs2 %d\n", tid, step->seq1[i].name, n_regs1, n_regs2);
	if (n_regs1 > 0)
	{
		assign_score (refFasta, &(step->seq1[i]), regs[0], n_regs1);
		qsort ((void *)regs[0], n_regs1, sizeof(mm_reg1_t), reg_compare);
//		printf ("\t%s score1 %d\n", step->seq1[i].name, regs[0][0].score);
	}
	if (n_regs2 > 0)
	{
		assign_score (refFasta, &(step->seq2[i]), regs[1], n_regs2);
		qsort ((void *)regs[1], n_regs2, sizeof(mm_reg1_t), reg_compare);
//		printf ("\t%s score2 %d\n", step->seq2[i].name, regs[1][0].score);
	}

	int mq1=-1, mq2=-1;
	if ((n_regs1 > 0) && (regs[0][0].score >= (step->seq1[i].l_seq - 4*4)))
	{
		int lseq2=0;
		if (step->seq2 != NULL) lseq2=step->seq2[i].l_seq;
		int regs2Ind = selectAlnRecs (regs[0], regs[1], n_regs1, n_regs2, step->seq1[i].l_seq, lseq2, &mq1, &mq2);
		if (regs2Ind != -1)
		{
//			if (regs[0][0].rid == 1)
//				printf ("1 - %s n_regs1 %d, n_regs2 %d - MQ1 %d MQ2 %d\n", step->seq1[i].name, n_regs1, n_regs2, regs[0][0].mapQ, regs[1][regs2Ind].mapQ);
			//vasu

			step->reg[2*i] = (mm_reg1_t*)malloc(sizeof(mm_reg1_t));
			memcpy(step->reg[2*i], regs[0], sizeof(mm_reg1_t));
			step->reg[2*i +1] = (mm_reg1_t*)malloc(sizeof(mm_reg1_t));
			memcpy(step->reg[2*i +1], &(regs[1][regs2Ind]), sizeof(mm_reg1_t));
			step->n_reg[2*i] = step->n_reg[2*i+1] = 1;

			fullRec fR;
			mm_reg1_t *r = step->reg[2*i];
//			r->fR = (fullRec *)malloc (sizeof(fullRec));
//			assert (r->fR != NULL);
//			printf ("call create_aeb-1 i %ld step->seq1 %p\n", i, step->seq1);
			create_aeb (&(step->seq1[i]), r, &fR);

			int32_t actRefPos = r->rs - r->qs;
			if (r->rev) actRefPos = r->rs - (step->seq1[i].l_seq - r->qe);
			if (actRefPos < 0) actRefPos = 0;
			int fileNo = r->rid*numMaxChunks + actRefPos/maxChunkSize;
			int oldSize = vfr[fileNo].m;
			kv_push_fr (vfr[fileNo], fR);
//			if (vfr[fileNo].m != oldSize)
//				printf ("%d: %d-->%d\n", fileNo, oldSize, vfr[fileNo].m);

			r = step->reg[2*i +1];
//			r->fR = (fullRec *)malloc (sizeof(fullRec));
//			assert (r->fR != NULL);
//			printf ("call create_aeb-2 i %ld\n", i);
			create_aeb (&(step->seq2[i]), r, &fR);

			actRefPos = r->rs - r->qs;
			if (r->rev) actRefPos = r->rs - (step->seq1[i].l_seq - r->qe);
			if (actRefPos < 0) actRefPos = 0;
			fileNo = r->rid*numMaxChunks + actRefPos/maxChunkSize;
			oldSize = vfr[fileNo].m;
			kv_push_fr (vfr[fileNo], fR);
//			if (vfr[fileNo].m != oldSize)
//				printf ("%d: %d-->%d\n", fileNo, oldSize, vfr[fileNo].m);
		}
		else
		{
//			printf ("Else1 %s (%d, %d)\n", step->seq1[i].name, step->n_reg[2*i], step->n_reg[2*i+1]);
			if (step->seq2 == NULL)
			{
				step->reg[2*i] = (mm_reg1_t*)malloc(sizeof(mm_reg1_t));
				memcpy(step->reg[2*i], regs[0], sizeof(mm_reg1_t));
				step->n_reg[2*i]=1;

				fullRec fR;
				mm_reg1_t *r = step->reg[2*i];
//				r->fR = (fullRec *)malloc (sizeof(fullRec));
//				assert (r->fR != NULL);
//				printf ("call create_aeb-3 i %ld\n", i);
				create_aeb (&(step->seq1[i]), r, &fR);

				int actRefPos = r->rs - r->qs;
				if (r->rev) actRefPos = r->rs - (step->seq1[i].l_seq - r->qe);
				if (actRefPos < 0) actRefPos = 0;
				int fileNo = r->rid*numMaxChunks + actRefPos/maxChunkSize;
				int oldSize = vfr[fileNo].m;
				kv_push_fr (vfr[fileNo], fR);
//				if (vfr[fileNo].m != oldSize)
//					printf ("%d: %d-->%d\n", fileNo, oldSize, vfr[fileNo].m);
			}
			else
			{//unmapped1
				int readSize = step->seq1[i].l_name + 3 + 2*step->seq1[i].l_seq + 1 + 4;
				if (step->seq2 != NULL)
					readSize += step->seq2[i].l_name + 3 + 2*step->seq2[i].l_seq + 1 + 4;
				if (step->p->u_n[tid] + readSize >= step->p->u_m[tid])
				{
					step->p->u_m[tid] = step->p->u_m[tid] ? step->p->u_m[tid] << 1 : 2*readSize;
					step->p->unmapped[tid] = (char *)realloc (step->p->unmapped[tid], step->p->u_m[tid]);
					assert (step->p->unmapped[tid] != NULL);
					bzero (step->p->unmapped[tid]+step->p->u_n[tid], step->p->u_m[tid]-step->p->u_n[tid]);
//					printf ("unmapped[%d]=%p size %lu \n", tid, step->p->unmapped[tid], step->p->u_m[tid]);
				}
				sprintf (&step->p->unmapped[tid][step->p->u_n[tid]], "@%s/1\n%s\n+\n%s\n", step->seq1[i].name, step->seq1[i].seq, step->seq1[i].qual);
				step->p->u_n[tid] += step->seq1[i].l_name + 3 + 2*step->seq1[i].l_seq + 1 + 4;
//				printf ("read %s read1Size %d\n", step->seq1[i].name, step->seq1[i].l_name + 3 + 2*step->seq1[i].l_seq + 1 + 4);
				if (step->seq2 != NULL)
				{
					sprintf (&step->p->unmapped[tid][step->p->u_n[tid]], "@%s/2\n%s\n+\n%s\n", step->seq2[i].name, step->seq2[i].seq, step->seq2[i].qual);
					step->p->u_n[tid] += step->seq2[i].l_name + 3 + 2*step->seq2[i].l_seq + 1 + 4;
				}

				step->n_reg[2*i] = step->n_reg[2*i+1] = 0;
			}
		}
	}
	else if ((n_regs2 > 0) && (regs[1][0].score >= step->seq2[i].l_seq - 4*4))
	{
		int regs1Ind = selectAlnRecs (regs[1], regs[0], n_regs2, n_regs1, step->seq2[i].l_seq, step->seq1[i].l_seq, &mq2, &mq1);
		if (regs1Ind != -1)
		{
//			if (regs[1][0].rid == 1)
//				printf ("2 - %s n_regs1 %d, n_regs2 %d - MQ1 %d MQ2 %d\n", step->seq1[i].name, n_regs1, n_regs2, mq1, mq2);
			step->reg[2*i] = (mm_reg1_t*)malloc(sizeof(mm_reg1_t));
			memcpy(step->reg[2*i], &(regs[0][regs1Ind]), sizeof(mm_reg1_t));
			step->reg[2*i +1] = (mm_reg1_t*)malloc(sizeof(mm_reg1_t));
			memcpy(step->reg[2*i +1], regs[1], sizeof(mm_reg1_t));
			step->n_reg[2*i] = step->n_reg[2*i+1] = 1;

			fullRec fR;
			mm_reg1_t *r = step->reg[2*i];
//			r->fR = (fullRec *)malloc (sizeof(fullRec));
//			assert (r->fR != NULL);
//			printf ("call create_aeb-4 i %ld\n", i);
			create_aeb (&(step->seq1[i]), r, &fR);

			int actRefPos = r->rs - r->qs;
			if (r->rev) actRefPos = r->rs - (step->seq1[i].l_seq - r->qe);
			if (actRefPos < 0) actRefPos = 0;
			int fileNo = r->rid*numMaxChunks + actRefPos/maxChunkSize;
			int oldSize = vfr[fileNo].m;
			kv_push_fr (vfr[fileNo], fR);
//			if (vfr[fileNo].m != oldSize)
//				printf ("%d: %d-->%d\n", fileNo, oldSize, vfr[fileNo].m);

			r = step->reg[2*i +1];
//			r->fR = (fullRec *)malloc (sizeof(fullRec));
//			assert (r->fR != NULL);
//			printf ("call create_aeb-5 i %ld\n", i);
			create_aeb (&(step->seq2[i]), r, &fR);

			actRefPos = r->rs - r->qs;
			if (r->rev) actRefPos = r->rs - (step->seq1[i].l_seq - r->qe);
			if (actRefPos < 0) actRefPos = 0;
			fileNo = r->rid*numMaxChunks + actRefPos/maxChunkSize;
			oldSize = vfr[fileNo].m;
			kv_push_fr (vfr[fileNo], fR);
//			if (vfr[fileNo].m != oldSize)
//				printf ("%d: %d-->%d\n", fileNo, oldSize, vfr[fileNo].m);
		}
		else
		{//unmapped2

			int readSize = step->seq1[i].l_name + 3 + 2*step->seq1[i].l_seq + 1 + 4;
			if (step->seq2 != NULL)
				readSize += step->seq2[i].l_name + 3 + 2*step->seq2[i].l_seq + 1 + 4;
			if (step->p->u_n[tid] + readSize >= step->p->u_m[tid])
			{
				step->p->u_m[tid] = step->p->u_m[tid] ? step->p->u_m[tid] << 1 : 2*readSize;
				step->p->unmapped[tid] = (char *)realloc (step->p->unmapped[tid], step->p->u_m[tid]);
				assert (step->p->unmapped[tid] != NULL);
				bzero (step->p->unmapped[tid]+step->p->u_n[tid], step->p->u_m[tid]-step->p->u_n[tid]);
//				printf ("unmapped[%d]=%p size %lu readLen %d\n", tid, step->p->unmapped[tid], step->p->u_m[tid], readSize);
			}
//			printf ("write *@%s/1*\n*%s*\n*+*\n*%s*\n*@%s/2*\n*%s*\n*+*\n*%s*\n(Len %d+%d) at pos %d,%d TotalSize %d\n", step->seq1[i].name, step->seq1[i].seq, step->seq1[i].qual, step->seq2[i].name, step->seq2[i].seq, step->seq2[i].qual, step->seq1[i].l_name + 3 + 2*step->seq1[i].l_seq + 1 + 4, step->seq2[i].l_name + 3 + 2*step->seq2[i].l_seq + 1 + 4, step->p->u_n[tid], step->p->u_n[tid]+step->seq1[i].l_name + 3 + 2*step->seq1[i].l_seq + 1 + 4, step->p->u_m[tid]);
			sprintf (&step->p->unmapped[tid][step->p->u_n[tid]], "@%s/1\n%s\n+\n%s\n", step->seq1[i].name, step->seq1[i].seq, step->seq1[i].qual);
			step->p->u_n[tid] += step->seq1[i].l_name + 3 + 2*step->seq1[i].l_seq + 1 + 4;
			if (step->seq2 != NULL)
			{
				sprintf (&step->p->unmapped[tid][step->p->u_n[tid]], "@%s/2\n%s\n+\n%s\n", step->seq2[i].name, step->seq2[i].seq, step->seq2[i].qual);
				step->p->u_n[tid] += step->seq2[i].l_name + 3 + 2*step->seq2[i].l_seq + 1 + 4;
			}

//			printf ("Else2 %s (%d, %d)\n", step->seq2[i].name, step->n_reg[2*i], step->n_reg[2*i+1]);
			step->n_reg[2*i] = step->n_reg[2*i+1] = 0;
		}
	}
	else //unmapped3 or complex read
	{

		int readSize = step->seq1[i].l_name + 3 + 2*step->seq1[i].l_seq + 1 + 4;
		if (step->seq2 != NULL)
			readSize += step->seq2[i].l_name + 3 + 2*step->seq2[i].l_seq + 1 + 4;
		if (step->p->u_n[tid] + readSize >= step->p->u_m[tid])
		{
			step->p->u_m[tid] = step->p->u_m[tid] ? step->p->u_m[tid] << 1 : 2*readSize;
			step->p->unmapped[tid] = (char *)realloc (step->p->unmapped[tid], step->p->u_m[tid]);
			assert (step->p->unmapped[tid] != NULL);
//			printf ("unmapped[%d]=%p size %lu \n", tid, step->p->unmapped[tid], step->p->u_m[tid]);
			bzero (step->p->unmapped[tid]+step->p->u_n[tid], step->p->u_m[tid]-step->p->u_n[tid]);
		}
		sprintf (&step->p->unmapped[tid][step->p->u_n[tid]], "@%s/1\n%s\n+\n%s\n", step->seq1[i].name, step->seq1[i].seq, step->seq1[i].qual);
		step->p->u_n[tid] += step->seq1[i].l_name + 3 + 2*step->seq1[i].l_seq + 1 + 4;
		if (step->seq2 != NULL)
		{
			sprintf (&step->p->unmapped[tid][step->p->u_n[tid]], "@%s/2\n%s\n+\n%s\n", step->seq2[i].name, step->seq2[i].seq, step->seq2[i].qual);
			step->p->u_n[tid] += step->seq2[i].l_name + 3 + 2*step->seq2[i].l_seq + 1 + 4;
		}

//		printf ("unmapped or complex read %s (%d, %d)\n", step->seq1[i].name, step->n_reg[2*i], step->n_reg[2*i+1]);
		step->n_reg[2*i] = step->n_reg[2*i+1] = 0;
	}
	
/*	if (i<=5)
		printf ("n_regs1 %d n_regs2 %d\n", n_regs1, n_regs2);*/
	if ((n_regs1 > 0) && ((step->seq2 == NULL) || (n_regs2 > 0))) {
/*		step->reg[2*i] = (mm_reg1_t*)malloc(n_regs1 * sizeof(mm_reg1_t));
		memcpy(step->reg[2*i], regs[0], n_regs1 * sizeof(mm_reg1_t));
		step->reg[2*i +1] = (mm_reg1_t*)malloc(n_regs2 * sizeof(mm_reg1_t));
		memcpy(step->reg[2*i +1], regs[1], n_regs2 * sizeof(mm_reg1_t));*/
	}
	else
	{
		step->n_reg[2*i] = step->n_reg[2*i+1] = 0;
/*		if (fpUnMapped[tid] == NULL)//write2
		{
			char fName[200];
			sprintf (fName, "/i3c/sparse/vxr162/minimap_pe/unmapped_%d_2.fq", tid);
			fpUnMapped[tid] = fopen (fName, "w");
		}
		assert (fpUnMapped[tid] != NULL);
		fprintf (fpUnMapped[tid], "@%s/1\n%s\n+\n%s\n", step->seq1[i].name, step->seq1[i].seq, step->seq1[i].qual);
		if (step->seq2 != NULL)
			fprintf (fpUnMapped[tid], "@%s/2\n%s\n+\n%s\n", step->seq2[i].name, step->seq2[i].seq, step->seq2[i].qual);*/
//		printf ("%d-%ld\n", tid, i);
	}
}
typedef struct read_thread
{
	pipeline_t *p;
	step_t *s;
	int tid;
}read_thread_t;
static void *map_worker(void *data)
{
	read_thread_t *rt = (read_thread_t *)data;
	pipeline_t *p = rt->p;
	step_t *s = rt->s;
	int tid = rt->tid;
	int i=0;
	step_t *in = (step_t*)calloc(1, sizeof(step_t));
	memcpy (in, s, sizeof(step_t));
	in->seq1 = s->seqs1[tid];
	in->seq2 = s->seqs2[tid];
	in->n_reg = s->n_regs[tid];
	in->reg = s->regs[tid];

	for (i=0; i<s->n_seq[tid]; i++)
	{
		worker_for(in, (long) i, tid);
	}
	free (in);
	return NULL;
}

static void *read_worker(void *data)
{
	read_thread_t *rt = (read_thread_t *)data;
	pipeline_t *p = rt->p;
	step_t *s = rt->s;
	int tid = rt->tid, i=0;
	
	double sTime = realtime();
	s->seqs1[tid] = bulk_read(p->fp1[tid], p->sizeToRead[tid], &s->n_seq[tid], &(p->sizeToRead[tid]), &(p->fastqBuffer1[tid]), &p->fastqSize1[tid]);
//	printf ("Tid %d FileOfs %d Size to read %d Time %lf sec nReads read %d\n", tid, bseq_tell (p->fp1[tid]), p->sizeToRead[tid], realtime()-sTime, s->n_seq[tid]);
/*	char fname[100];
	sprintf (fname, "%s/thr%d.fastq", p->oprefix, tid);
	FILE *fout = fopen (fname, "a");
	assert (fout != NULL);
	fwrite (s->seqs1[tid], 1, readSize, fout);
	fclose (fout);
	free (s->seqs1[tid]);
	s->seqs1[tid] = NULL;*/

	if (tid == 0)
	{
		s->buf1 = (mm_tbuf_t**)calloc(p->n_threads, sizeof(mm_tbuf_t*));
//		printf ("CALLOC %p\n", s->buf1);
		for (i = 0; i < p->n_threads; ++i)
			s->buf1[i] = mm_tbuf_init();
	}
	if ((s->seqs1[tid] != NULL) && (p->fp2[tid] != NULL))
	{
//		s->seqs2[tid] = bseq_read2(p->fp2[tid], s->n_seq[tid]);
		int n;
		s->seqs2[tid] = bulk_read(p->fp2[tid], p->sizeToRead2[tid], &n, &(p->sizeToRead[tid]), &(p->fastqBuffer2[tid]), &p->fastqSize2[tid]);
		assert ((s->seqs2[tid] != NULL) && (n == s->n_seq[tid]));
		if (tid == 0)
		{
			s->buf2 = (mm_tbuf_t**)calloc(p->n_threads, sizeof(mm_tbuf_t*));
			for (i = 0; i < p->n_threads; ++i)
				s->buf2[i] = mm_tbuf_init();
		}
	}
//	if (s->seqs1[tid] != NULL)
//		printf ("After Tid %d FileOfs %lu numSeq %d name %s\n", tid, bseq_tell (p->fp1[tid]), s->n_seq[tid], s->seqs1[tid][0].name);
	pthread_exit(0);
}
static void *worker_pipeline(void *shared, int step, void *in)
{
	int i, j, k;
    pipeline_t *p = (pipeline_t*)shared;
    int numOutChunks = p->numOutChunks;
    int maxChunkSize = p->maxChunkSize;
    int numMaxChunks = p->numMaxChunks;
    FILE **aebFp = p->aebFp;
//    printf ("worker_pipeline p->refFasta %p p->mi %p numOutChunks %d, maxChunkSize %d, numMaxChunks %d\n", p->refFasta, p->mi, numOutChunks, maxChunkSize, numMaxChunks);
    double sTime = realtime();
//	printf ("\tworker_pipeline Start step %d at %lf sec\n", step, sTime);
    double sTime_cpu = cputime();
    if (step == 0) { // step 0: read sequences
        step_t *s;
        s = (step_t*)calloc(1, sizeof(step_t));
	s->seqs1 = (bseq1_t **)calloc (p->n_threads, sizeof(bseq1_t *));
	s->seqs2 = (bseq1_t **)calloc (p->n_threads, sizeof(bseq1_t *));
	s->n_seq = (int *)calloc (p->n_threads, sizeof(int));

	s->n_regs = (int **)calloc (p->n_threads, sizeof(int *));
	s->regs = (mm_reg1_t***) calloc(p->n_threads, sizeof(mm_reg1_t**));
	s->p = p;

	int tid=0;
	read_thread_t *rt = (read_thread_t *)malloc (p->n_threads * sizeof(read_thread_t));
	assert (rt != NULL);
	pthread_t *thr = (pthread_t *)malloc (p->n_threads * sizeof (pthread_t));
	assert (thr != NULL);
	for (tid = 0; tid<p->n_threads; tid++)
	{
		rt[tid].tid = tid;
		rt[tid].p = p;
		rt[tid].s = s;
		pthread_create (&thr[tid], 0, read_worker, &rt[tid]);
	}
	for (tid = 0; tid<p->n_threads; tid++)
		pthread_join (thr[tid], 0);
	free (thr);
	free (rt);
//	printf ("\t1-worker_pipeline step %d(fileRead) real %.3lf sec cpu %.3lf sec\n",  step, realtime()-sTime, cputime()-sTime_cpu);

	int rec_found=0;
	for (tid = 0; tid<p->n_threads; tid++)
	{
		if (s->seqs1[tid]) {
			rec_found=1;
			if (p->fp2[tid] != NULL)
			{
				int n_seq=0;
//				sTime = realtime();
				for (i = 0; i < s->n_seq[tid]; ++i)
				{
//					printf ("tid %d read name %s nameLen %d seq %s seqLen %d Qual %s\n", tid, s->seqs1[tid][i].name, s->seqs1[tid][i].l_name, s->seqs1[tid][i].seq, s->seqs1[tid][i].l_seq, s->seqs1[tid][i].qual);
					s->seqs1[tid][i].rid = p->n_processed++;
					s->seqs2[tid][i].rid = p->n_processed++;
				}
			}
			else
			{
				for (i = 0; i < s->n_seq[tid]; ++i)
					s->seqs1[tid][i].rid = p->n_processed++;
				s->seqs2[tid] = NULL;
			}
			s->n_regs[tid] = (int *)calloc(s->n_seq[tid]*2, sizeof(int));
			s->regs[tid] = (mm_reg1_t**)calloc(s->n_seq[tid]*2, sizeof(mm_reg1_t*));

		} else {
		}
	}
	printf ("\t1-worker_pipeline step %d real %.3lf sec cpu %.3lf sec recFound %d\n",  step, realtime()-sTime, cputime()-sTime_cpu, rec_found);
	if (rec_found)
		return s;
	else
	{
		for (i = 0; i < p->n_threads; ++i) mm_tbuf_destroy(s->buf1[i]);//to free
//		printf ("FREE2 %p\n", s->buf1);
		free(s->buf1);
		if (s->seqs2[0] != NULL)
		{
			for (i = 0; i < p->n_threads; ++i) mm_tbuf_destroy(s->buf2[i]);
			free(s->buf2);
		}
		free (s->seqs1);
		free (s->seqs2);
		free (s->n_regs);
		free (s->regs);
		free (s->n_seq);
		free(s);
	}
    } else if (step == 1) { // step 1: map
/*		for (i=0; i< p->n_threads; i++)
		{
			((step_t*)in)->seq1 = ((step_t*)in)->seqs1[i];
			((step_t*)in)->seq2 = ((step_t*)in)->seqs2[i];
			((step_t*)in)->n_reg = ((step_t*)in)->n_regs[i];
			((step_t*)in)->reg = ((step_t*)in)->regs[i];
//			printf ("tid %d step1 %p step2 %p\n", i, ((step_t*)in)->seq1, ((step_t*)in)->seq2);
			kt_for(p->n_threads, worker_for, in, ((step_t*)in)->n_seq[i]);
		}*/
		read_thread_t *rt = (read_thread_t *)malloc (p->n_threads * sizeof(read_thread_t));
		assert (rt != NULL);
		pthread_t *thr = (pthread_t *)malloc (p->n_threads * sizeof (pthread_t));
		assert (thr != NULL);
		for (i = 0; i<p->n_threads; i++)
		{
			rt[i].tid = i;
			rt[i].p = p;
			rt[i].s = (step_t *)in;
			pthread_create (&thr[i], 0, map_worker, &rt[i]);
		}
		for (i = 0; i<p->n_threads; i++)
			pthread_join (thr[i], 0);
		free (thr);
		free (rt);

    		printf ("\t2-worker_pipeline step %d real %.3lf sec cpu %.3lf sec\n",  step, realtime()-sTime, cputime()-sTime_cpu);
		return in;
    } else if (step == 2) { // step 2: output
        step_t *s = (step_t*)in;
//	printf ("STEP2\n");

	const mm_idx_t *mi = p->mi;
	for (i = 0; i < p->n_threads; ++i) mm_tbuf_destroy(s->buf1[i]);//to free
//	printf ("FREE %p\n", s->buf1);
	free(s->buf1);
	if (s->seqs2[0] != NULL)
	{
		for (i = 0; i < p->n_threads; ++i) mm_tbuf_destroy(s->buf2[i]);
		free(s->buf2);
	}
	for (i=0; i<p->n_threads; i++)
	{
		for (j=0; j<p->mi->n; j++)
		{
			for (k=0; k<numMaxChunks; k++)
			{
				int vecNo = i * p->mi->n * numMaxChunks + j * numMaxChunks + k;
				int fileNo = j * numMaxChunks + k;
				if (p->vfr[vecNo].n > 0) 
				{
					if (aebFp[fileNo] == NULL)
					{
						char fName[500];
						sprintf (fName, "%s/C%d_%d.aeb", p->oprefix, j, k);
						aebFp[fileNo] = fopen (fName, "w");
						assert (aebFp[fileNo] != NULL);
					}
//					printf ("vecNo %d FileNo %d - T%d C%d F%d : n:%lu a:%p\n", vecNo, fileNo, i, j, k, p->vfr[vecNo].n, p->vfr[vecNo].a);
					fwrite (p->vfr[vecNo].a, sizeof (fullRec), p->vfr[vecNo].n, aebFp[fileNo]);//aeb write
					free (p->vfr[vecNo].a);
					p->vfr[vecNo].a = 0;
					p->vfr[vecNo].n = 0;
					p->vfr[vecNo].m = 0;
				}
			}
		}
		if (p->u_n[i] > 0)
		{
			if (aebFp[p->mi->n*numMaxChunks] == NULL)//write2
			{
				char fName[500];
				sprintf (fName, "%s/unmapped.fq", p->oprefix);
				aebFp[p->mi->n*numMaxChunks] = fopen (fName, "w");
			}
			assert (aebFp[p->mi->n*numMaxChunks] != NULL);
			fwrite (p->unmapped[i], 1, p->u_n[i], aebFp[p->mi->n*numMaxChunks]);
//			printf ("unmapped %d \n%s\n", i, p->unmapped[i]);
			p->u_n[i] = 0;
		}
	}

	int tid=0;
	for (tid=0; tid<p->n_threads; tid++)
	{//tid loop
		for (i = 0; i < s->n_seq[tid]; ++i) {
			bseq1_t *t = &s->seqs1[tid][i];
			int found=0;
	                if ((t->rid % 1000000) == 0)
	                        fprintf (stderr, "%d\t%p\n", t->rid, t->name);
	/*		if ((s->n_reg[2*i]== 0) && (s->n_reg[2*i +1] == 0))
			{
				if (aebFp[p->mi->n*numMaxChunks] == NULL)//write2
				{
					char fName[500];
					sprintf (fName, "%s/unmapped.fq", p->oprefix);
					aebFp[p->mi->n*numMaxChunks] = fopen (fName, "w");
				}
				assert (aebFp[p->mi->n*numMaxChunks] != NULL);
				fprintf (aebFp[p->mi->n*numMaxChunks], "@%s/1\n%s\n+\n%s\n", s->seq1[i].name, s->seq1[i].seq, s->seq1[i].qual);
				if (s->seq2 != NULL)
					fprintf (aebFp[p->mi->n*numMaxChunks], "@%s/2\n%s\n+\n%s\n", s->seq2[i].name, s->seq2[i].seq, s->seq2[i].qual);
			}*/
	/*		fullRec fR;
			for (j = 0; j < s->n_reg[2*i]; ++j) {
				mm_reg1_t *r = &s->reg[2*i][j];
	
				int32_t actRefPos=r->rs - r->qs;
				if (r->rev) actRefPos = r->rs - (t->l_seq - r->qe);
				if (actRefPos < 0) actRefPos = 0;
				int fileNo = r->rid*numMaxChunks + actRefPos/maxChunkSize;
				if (aebFp[fileNo] == NULL)
				{
					char fName[500];
					sprintf (fName, "%s/C%d_%d.aeb", p->oprefix, r->rid, actRefPos/maxChunkSize);
					aebFp[fileNo] = fopen (fName, "w");
					assert (aebFp[fileNo] != NULL);
				}
				fwrite (r->fR, sizeof (fullRec), 1, aebFp[fileNo]);
				free (r->fR);
				r->fR=NULL;
	//			if (r->len < p->opt->min_match) continue;//vasu
	
	//				printf("%d\t%s\t%d\t%s\t%s\t%d\t%d\t%c\t", t->rid, t->name, t->l_seq, t->seq, t->qual, r->qs, r->qe, "+-"[r->rev]);
	//				printf("%d", r->rid);
	//				printf("\t%d\t%d\t%d\t%d\t%d\t255\tcm:i:%d\tscore %d\n", mi->len[r->rid], r->rs, r->re, r->len,
	//						r->re - r->rs > r->qe - r->qs? r->re - r->rs : r->qe - r->qs, r->cnt, r->score);
			}*/
			if (s->n_regs[tid][2*i] > 0)
				free(s->regs[tid][2*i]);
			if (s->seqs2[tid] != NULL)
	      		{
				t = &s->seqs2[tid][i];
	/*			for (j = 0; j < s->n_reg[2*i+1]; ++j) {
					mm_reg1_t *r = &s->reg[2*i+1][j];
					int32_t actRefPos=r->rs - r->qs;
					if (r->rev) actRefPos = r->rs - (t->l_seq - r->qe);
					if (actRefPos < 0) actRefPos=0;
					int fileNo = r->rid*numMaxChunks + actRefPos/maxChunkSize;
					if (aebFp[fileNo] == NULL)
					{
						char fName[500];
						sprintf (fName, "%s/C%d_%d.aeb", p->oprefix, r->rid, actRefPos/maxChunkSize);
						aebFp[fileNo] = fopen (fName, "w");
						assert (aebFp[fileNo] != NULL);
					}
					fwrite (r->fR, sizeof (fullRec), 1, aebFp[fileNo]);
					free (r->fR);
					r->fR=NULL;
	//				if (r->len < p->opt->min_match) continue;//vasu
	
	//					printf("%d\t%s\t%d\t%s\t%s\t%d\t%d\t%c\t", t->rid, t->name, t->l_seq, t->seq, t->qual, r->qs, r->qe, "+-"[r->rev]);
	//					printf("%d", r->rid);
	//					printf("\t%d\t%d\t%d\t%d\t%d\t255\tcm:i:%d\tscore %d\n", mi->len[r->rid], r->rs, r->re, r->len,
	//							r->re - r->rs > r->qe - r->qs? r->re - r->rs : r->qe - r->qs, r->cnt, r->score);
				}*/
				if (s->n_regs[tid][2*i +1] > 0)
					free(s->regs[tid][2*i +1]);
			}
	//		printf ("s->seq1.name %p %p\n", s->seq1, s->seq1[i]);
			if (i==0)
			{
	//			free(s->seqs1[tid][i].name-1); //frees block of fastq data (all n_seq[tid])
	//			printf ("free %p\n", s->seqs1[tid][i].name-1);
				if (s->seqs2[tid])
				{
	//				free(s->seqs2[tid][i].name-1); //frees block of fastq data (all n_seq[tid])
	//				printf ("free %p\n", s->seqs2[tid][i].name-1);
				}
			}
	/*		free(s->seqs1[tid][i].seq); 
			free(s->seqs1[tid][i].name); 
			free (s->seqs1[tid][i].qual);
			if (s->seqs2[tid])
			{
				free(s->seqs2[tid][i].seq); free(s->seqs2[tid][i].name); free (s->seqs2[tid][i].qual);
			}*/
		}
//		printf ("free seqs1 %p \n", s->seqs1[tid]);
		free(s->seqs1[tid]); free (s->n_regs[tid]); free (s->regs[tid]);
		if (s->seqs2[tid]) 
		{
//			printf ("free seqs2 %p \n", s->seqs2[tid]);
			free(s->seqs2[tid]);
		}
	}//tid loop

//		printf ("step2 out\n");
	free(s->regs); free(s->n_regs); free(s->seqs1); 
	free (s->n_seq);
	free(s->seqs2);
	free(s);
    }
    printf ("\t3-worker_pipeline step %d real %.3lf sec cpu %.3lf sec\n",  step, realtime()-sTime, cputime()-sTime_cpu);
    return 0;
}

int mm_map_file(int numTasks, int rank, int numInChunks, int numOutChunks, int maxChunkSize, int numMaxChunks, fileIndex *fI, FILE **aebFp, char *oprefix, char **refFasta, const mm_idx_t *idx, const char *fn1, const char *fn2, const mm_mapopt_t *opt, int n_threads, int tbatch_size)
{
//	fI[1]=fI[numInChunks-1];
	pipeline_t pl;
	int i=0, step=0;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.fp1 = (bseq_file_t **) calloc (n_threads, sizeof(bseq_file_t *));
	pl.fp2 = (bseq_file_t **) calloc (n_threads, sizeof(bseq_file_t *));
	pl.sizeToRead = (int *)calloc (n_threads, sizeof(int));
	pl.sizeToRead2 = (int *)calloc (n_threads, sizeof(int));
	pl.fastqBuffer1 = (char **)calloc (n_threads, sizeof(char *));
	pl.fastqBuffer2 = (char **)calloc (n_threads, sizeof(char *));
	pl.fastqSize1 = (size_t *)calloc (n_threads, sizeof(size_t));
	pl.fastqSize2 = (size_t *)calloc (n_threads, sizeof(size_t));
	assert ((pl.fp1 != NULL) && (pl.fp2 != NULL) && (pl.sizeToRead != NULL) && (pl.sizeToRead2 != NULL) && (pl.fastqBuffer1 != NULL) && (pl.fastqBuffer2 != NULL) && (pl.fastqSize1 != NULL) && (pl.fastqSize2 != NULL));

	pl.refFasta = refFasta;
	pl.oprefix = oprefix;

//	kvec_t(fullRec) tmp;

	vec_fullRec *vfr = (vec_fullRec *)malloc(n_threads * numMaxChunks * idx->n * sizeof(vec_fullRec));
	assert (vfr != NULL);
//	vfr++;
	for (i=0;i<n_threads*numMaxChunks * idx->n; i++)
		kv_init(vfr[i]);

	pl.numOutChunks = numOutChunks;
	pl.maxChunkSize = maxChunkSize;
	pl.numMaxChunks = numMaxChunks;
	pl.aebFp = aebFp;
	pl.vfr = vfr;
	pl.unmapped = (char **)calloc (n_threads, sizeof(char *));
	pl.u_m = (size_t *)calloc (n_threads, sizeof (size_t));
	pl.u_n = (size_t *)calloc (n_threads, sizeof (size_t));
	assert ((pl.unmapped != NULL) && (pl.u_m != NULL) && (pl.u_n != NULL));
	pl.opt = opt, pl.mi = idx;
	pl.n_threads = n_threads, pl.batch_size = tbatch_size;
//	printf ("n_threads %d\n", pl.n_threads);
	
//	for (;i<n_threads; i++)
		fpUnMapped=NULL;
	double rTime1=realtime();
	for (i=0; i<n_threads; i++)
	{
	}
	for (step=rank*n_threads; step <numInChunks-1; step+=numTasks*n_threads)
	{
		printf ("[%d]STEP_IDX %d/%d\n", rank, step, numInChunks);
		for (i=0; i<n_threads; i++)
		{
			pl.fp1[i] = bseq_open(fn1);//vasu
			if (pl.fp1 == 0) return -1;
			if (fn2 != NULL)
			{
				pl.fp2[i] = bseq_open(fn2);
				if (pl.fp2[i] == 0) return -1;
			}

			pl.sizeToRead[i] = fI[step+i+1].offset-fI[step+i].offset;
			pl.sizeToRead2[i] = fI[step+i+1].offset2-fI[step+i].offset2;
			bseq_seek (pl.fp1[i], fI[step+i].offset);
	//	printf ("fn1 %s fn2 %s\n", fn1, fn2);
			if (fn2 != NULL)
			{
				bseq_seek (pl.fp2[i], fI[step+i].offset2);
			}
		}
	
		kt_pipeline(n_threads == 1? 1 : 1, worker_pipeline, &pl, 3);

		for (i=0; i<n_threads; i++)
		{
			bseq_close(pl.fp1[i]);
			if (pl.fp2[i] != NULL) bseq_close(pl.fp2[i]);
		}
	}

	printf ("kt_pipeline time %.3lf sec\n", realtime()-rTime1);
//	for (i=0;i<n_threads; i++)
//	{
		if (fpUnMapped != NULL)
		{
			fclose (fpUnMapped);
			fpUnMapped = NULL;
		}
//	}
	for (i=0;i<n_threads*numMaxChunks * idx->n; i++)
	{
		if (vfr[i].a) {
			free (vfr[i].a);
			vfr[i].a = 0;
			vfr[i].m = 0;
			vfr[i].n = 0;
		}
	}
	for (i=0; i<n_threads; i++)
	{
		if (pl.unmapped[i] != NULL) free (pl.unmapped[i]);
		if (pl.fastqBuffer1[i]  != NULL) 
		{
			free (pl.fastqBuffer1[i]);
		}
		if (pl.fastqBuffer2[i]  != NULL) 
		{
			free (pl.fastqBuffer2[i]);
		}
//		printf ("Free unmapped[%d]=%p size %lu \n", i, pl.unmapped[i], pl.u_m[i]);
	}
	free (pl.unmapped);
	free (pl.u_m);
	free (pl.u_n);
	free (pl.fp1);
	free (pl.fp2);
	free (pl.sizeToRead);
	free (pl.sizeToRead2);
	free (pl.fastqBuffer1);
	free (pl.fastqBuffer2);
	free (pl.fastqSize1);
	free (pl.fastqSize2);
	free (vfr);
	return 0;
}
