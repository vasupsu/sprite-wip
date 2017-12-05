#include <zlib.h>
#ifdef USE_MPI
#include <mpi.h>
#endif
#include <pthread.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "bwa.h"
#include "bwamem.h"
#include "kvec.h"
#include "utils.h"
#include "kseq.h"
#include "utils.h"
KSEQ_DECLARE(gzFile)

extern unsigned char nst_nt4_table[256];

void *kopen(const char *fn, int *_fd);
int kclose(void *a);
char *bwa_pg;
struct stat buf;

typedef struct parsnip_data
{
        fullRec *fR;
        long len;
        int contig;
        int tid;
        long *nextOfs;
        long *remLen;
}parsnip_data;

typedef struct
{
   FILE *fp_fullMatch;
   FILE *fp_others;
   char *fullFileName;
   char *otherFileName;
   long curFullCount;
   long curOtherCount;
   long outputChunkSize;
   pthread_mutex_t lock;
   fullRec *Fr;    
   otherRec *Or;
}contigInfo;

extern void parsnip_init(int rank, int seqLen, char *fasta_file, char *vcf_prefix, int numContigs, int contig, int nt);
extern void *parsnip_aeb(void *data);
extern void parsnip_aib(otherRec *Or, long numOr, int contig);
extern void *parsnip_write(void *data);

struct timeval sTime, eTime;
long compE=0, overE=0, idxE=0, readE=0, writeE=0;
static void update_a(mem_opt_t *opt, const mem_opt_t *opt0)
{
	if (opt0->a) { // matching score is changed
		if (!opt0->b) opt->b *= opt->a;
		if (!opt0->T) opt->T *= opt->a;
		if (!opt0->o_del) opt->o_del *= opt->a;
		if (!opt0->e_del) opt->e_del *= opt->a;
		if (!opt0->o_ins) opt->o_ins *= opt->a;
		if (!opt0->e_ins) opt->e_ins *= opt->a;
		if (!opt0->zdrop) opt->zdrop *= opt->a;
		if (!opt0->pen_clip5) opt->pen_clip5 *= opt->a;
		if (!opt0->pen_clip3) opt->pen_clip3 *= opt->a;
		if (!opt0->pen_unpaired) opt->pen_unpaired *= opt->a;
	}
}

void bwa_print_sam_hdr2(const bntseq_t *bns, const char *rg_line, FILE *fp)
{
        int i;
        extern char *bwa_pg;
        for (i = 0; i < bns->n_seqs; ++i)
                fprintf(fp, "@SQ\tSN:%s\tLN:%d\n", bns->anns[i].name, bns->anns[i].len);
        if (rg_line) fprintf(fp, "%s\n", rg_line);
        fprintf(fp, "%s\n", bwa_pg);
}
#ifdef USE_MPI
void createAebType (MPI_Datatype *type)
{
	const int nitems=4;
        int blocklengths[4] = {1,26,1,1};
        MPI_Datatype types[4] = {MPI_UNSIGNED, MPI_BYTE, MPI_BYTE, MPI_BYTE};
        MPI_Aint offsets[4];
        offsets[0] = offsetof(fullRec, pos);
        offsets[1] = offsetof(fullRec, seq);
        offsets[2] = offsetof(fullRec, qual);
        offsets[3] = offsetof(fullRec, matchLen);

        MPI_Type_create_struct(nitems, blocklengths, offsets, types, type);
}
void createAibType (MPI_Datatype *type)
{
	const int nitems=5;
	int blocklengths[5] = {1,26,1,1,5};
	MPI_Datatype types[5] = {MPI_UNSIGNED, MPI_BYTE, MPI_BYTE, MPI_BYTE, MPI_UNSIGNED};
	MPI_Aint offsets[5];
        offsets[0] = offsetof(otherRec, pos);
        offsets[1] = offsetof(otherRec, seq);
        offsets[2] = offsetof(otherRec, qual);
        offsets[3] = offsetof(otherRec, n_cigar);
        offsets[4] = offsetof(otherRec, cigar);
	
	MPI_Type_create_struct(nitems, blocklengths, offsets, types, type);
}
#endif
static int cmpFullRec(const void *p1, const void *p2)
{
	fullRec *f1=(fullRec *)p1;
	fullRec *f2=(fullRec *)p2;
           return f1->pos-f2->pos;
}
static int cmpOtherRec(const void *p1, const void *p2)
{
	otherRec *f1=(otherRec *)p1;
	otherRec *f2=(otherRec *)p2;
           return f1->pos-f2->pos;
}
int main(int argc, char *argv[])
{
	mem_opt_t *opt, opt0;
	int fd, fd2, i, c, n, copy_comment = 0, j;
	gzFile fp, fp2 = 0;
	kseq_t *ks, *ks2 = 0;
	bseq1_t *seqs;
	bwaidx_t *idx;
	char *p, *rg_line = 0;
	const char *mode = 0;
	void *ko = 0, *ko2 = 0;
	int numtasks=1, rank=0;
	int64_t n_processed = 0;
	mem_pestat_t pes[4], *pes0 = 0;
	char *out_prefix=NULL;
#ifdef USE_MPI
	MPI_Init(&argc,&argv);
        MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif

	memset(pes, 0, 4 * sizeof(mem_pestat_t));
	for (i = 0; i < 4; ++i) pes[i].failed = 1;

	opt = mem_opt_init();
	memset(&opt0, 0, sizeof(mem_opt_t));
	while ((c = getopt(argc, argv, "epaFMCSPHYk:c:v:r:t:R:A:B:O:E:U:w:L:d:T:Q:D:m:I:N:W:x:G:h:o:")) >= 0) {
		if (c == 'k') opt->min_seed_len = atoi(optarg), opt0.min_seed_len = 1;
		else if (c == 'x') mode = optarg;
		else if (c == 'w') opt->w = atoi(optarg), opt0.w = 1;
		else if (c == 'A') opt->a = atoi(optarg), opt0.a = 1;
		else if (c == 'B') opt->b = atoi(optarg), opt0.b = 1;
		else if (c == 'T') opt->T = atoi(optarg), opt0.T = 1;
		else if (c == 'U') opt->pen_unpaired = atoi(optarg), opt0.pen_unpaired = 1;
		else if (c == 't') opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1;
		else if (c == 'P') opt->flag |= MEM_F_NOPAIRING;
		else if (c == 'a') opt->flag |= MEM_F_ALL;
		else if (c == 'p') opt->flag |= MEM_F_PE;
		else if (c == 'M') opt->flag |= MEM_F_NO_MULTI;
		else if (c == 'S') opt->flag |= MEM_F_NO_RESCUE;
		else if (c == 'e') opt->flag |= MEM_F_SELF_OVLP;
		else if (c == 'F') opt->flag |= MEM_F_ALN_REG;
		else if (c == 'Y') opt->flag |= MEM_F_SOFTCLIP;
		else if (c == 'c') opt->max_occ = atoi(optarg), opt0.max_occ = 1;
		else if (c == 'd') opt->zdrop = atoi(optarg), opt0.zdrop = 1;
		else if (c == 'v') bwa_verbose = atoi(optarg);
		else if (c == 'r') opt->split_factor = atof(optarg), opt0.split_factor = 1.;
		else if (c == 'D') opt->drop_ratio = atof(optarg), opt0.drop_ratio = 1.;
		else if (c == 'm') opt->max_matesw = atoi(optarg), opt0.max_matesw = 1;
//		else if (c == 'h') opt->max_hits = atoi(optarg), opt0.max_hits = 1;
		else if (c == 's') opt->split_width = atoi(optarg), opt0.split_width = 1;
		else if (c == 'G') opt->max_chain_gap = atoi(optarg), opt0.max_chain_gap = 1;
		else if (c == 'N') opt->max_chain_extend = atoi(optarg), opt0.max_chain_extend = 1;
		else if (c == 'W') opt->min_chain_weight = atoi(optarg), opt0.min_chain_weight = 1;
		else if (c == 'C') copy_comment = 1;
		else if (c == 'o') out_prefix=optarg;
		else if (c == 'Q') {
			opt0.mapQ_coef_len = 1;
			opt->mapQ_coef_len = atoi(optarg);
			opt->mapQ_coef_fac = opt->mapQ_coef_len > 0? log(opt->mapQ_coef_len) : 0;
		} else if (c == 'O') {
			opt0.o_del = opt0.o_ins = 1;
			opt->o_del = opt->o_ins = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->o_ins = strtol(p+1, &p, 10);
		} else if (c == 'E') {
			opt0.e_del = opt0.e_ins = 1;
			opt->e_del = opt->e_ins = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->e_ins = strtol(p+1, &p, 10);
		} else if (c == 'L') {
			opt0.pen_clip5 = opt0.pen_clip3 = 1;
			opt->pen_clip5 = opt->pen_clip3 = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->pen_clip3 = strtol(p+1, &p, 10);
		} else if (c == 'R') {
			if ((rg_line = bwa_set_rg(optarg)) == 0) return 1; // FIXME: memory leak
		} else if (c == 'I') { // specify the insert size distribution
			pes0 = pes;
			pes[1].failed = 0;
			pes[1].avg = strtod(optarg, &p);
			pes[1].std = pes[1].avg * .1;
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].std = strtod(p+1, &p);
			pes[1].high = (int)(pes[1].avg + 4. * pes[1].std + .499);
			pes[1].low  = (int)(pes[1].avg - 4. * pes[1].std + .499);
			if (pes[1].low < 1) pes[1].low = 1;
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].high = (int)(strtod(p+1, &p) + .499);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].low  = (int)(strtod(p+1, &p) + .499);
			if (bwa_verbose >= 3)
				fprintf(stderr, "[M::%s] mean insert size: %.3f, stddev: %.3f, max: %d, min: %d\n",
						__func__, pes[1].avg, pes[1].std, pes[1].high, pes[1].low);
		}
		else return 1;
	}
	if (opt->n_threads < 1) opt->n_threads = 1;
	if (optind + 1 >= argc || optind + 4 < argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq] [output prefix. Default: in1.fq]\n\n");
		fprintf(stderr, "Algorithm options:\n\n");
		fprintf(stderr, "       -t INT        number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "       -k INT        minimum seed length [%d]\n", opt->min_seed_len);
		fprintf(stderr, "       -w INT        band width for banded alignment [%d]\n", opt->w);
		fprintf(stderr, "       -d INT        off-diagonal X-dropoff [%d]\n", opt->zdrop);
		fprintf(stderr, "       -r FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [%g]\n", opt->split_factor);
//		fprintf(stderr, "       -s INT        look for internal seeds inside a seed with less than INT occ [%d]\n", opt->split_width);
		fprintf(stderr, "       -c INT        skip seeds with more than INT occurrences [%d]\n", opt->max_occ);
		fprintf(stderr, "       -D FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [%.2f]\n", opt->drop_ratio);
		fprintf(stderr, "       -W INT        discard a chain if seeded bases shorter than INT [0]\n");
		fprintf(stderr, "       -m INT        perform at most INT rounds of mate rescues for each read [%d]\n", opt->max_matesw);
		fprintf(stderr, "       -S            skip mate rescue\n");
		fprintf(stderr, "       -P            skip pairing; mate rescue performed unless -S also in use\n");
		fprintf(stderr, "       -e            discard full-length exact matches\n");
		fprintf(stderr, "       -A INT        score for a sequence match, which scales options -TdBOELU unless overridden [%d]\n", opt->a);
		fprintf(stderr, "       -B INT        penalty for a mismatch [%d]\n", opt->b);
		fprintf(stderr, "       -O INT[,INT]  gap open penalties for deletions and insertions [%d,%d]\n", opt->o_del, opt->o_ins);
		fprintf(stderr, "       -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [%d,%d]\n", opt->e_del, opt->e_ins);
		fprintf(stderr, "       -L INT[,INT]  penalty for 5'- and 3'-end clipping [%d,%d]\n", opt->pen_clip5, opt->pen_clip3);
		fprintf(stderr, "       -U INT        penalty for an unpaired read pair [%d]\n", opt->pen_unpaired);
		fprintf(stderr, "       -x STR        read type. Setting -x changes multiple parameters unless overriden [null]\n");
		fprintf(stderr, "                     pacbio: -k17 -W40 -r10 -A2 -B5 -O2 -E1 -L0\n");
		fprintf(stderr, "                     pbread: -k13 -W40 -c1000 -r10 -A2 -B5 -O2 -E1 -N25 -FeaD.001\n");
		fprintf(stderr, "\nInput/output options:\n\n");
		fprintf(stderr, "       -p            first query file consists of interleaved paired-end sequences\n");
		fprintf(stderr, "       -R STR        read group header line such as '@RG\\tID:foo\\tSM:bar' [null]\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "       -v INT        verbose level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", bwa_verbose);
		fprintf(stderr, "       -T INT        minimum score to output [%d]\n", opt->T);
//		fprintf(stderr, "       -h INT        if there are <INT hits with score >80%% of the max score, output all in XA [%d]\n", opt->max_hits);
		fprintf(stderr, "       -a            output all alignments for SE or unpaired PE\n");
		fprintf(stderr, "       -o STR        output file prefix (Default: in1.fq)\n");
		fprintf(stderr, "       -C            append FASTA/FASTQ comment to SAM output\n");
		fprintf(stderr, "       -Y            use soft clipping for supplementary alignments\n");
		fprintf(stderr, "       -M            mark shorter split hits as secondary\n\n");
		fprintf(stderr, "       -I FLOAT[,FLOAT[,INT[,INT]]]\n");
		fprintf(stderr, "                     specify the mean, standard deviation (10%% of the mean if absent), max\n");
		fprintf(stderr, "                     (4 sigma from the mean if absent) and min of the insert size distribution.\n");
		fprintf(stderr, "                     FR orientation only. [inferred]\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "Note: Please read the man page for detailed description of the command line and options.\n");
		fprintf(stderr, "\n");
		free(opt);
		return 1;
	}

	if (mode) {
		if (strcmp(mode, "pacbio") == 0 || strcmp(mode, "pbref") == 0 || strcmp(mode, "pbread1") == 0 || strcmp(mode, "pbread") == 0) {
			if (!opt0.a) opt->a = 2, opt0.a = 1;
			update_a(opt, &opt0);
			if (!opt0.o_del) opt->o_del = 2;
			if (!opt0.e_del) opt->e_del = 1;
			if (!opt0.o_ins) opt->o_ins = 2;
			if (!opt0.e_ins) opt->e_ins = 1;
			if (!opt0.b) opt->b = 5;
			if (opt0.split_factor == 0.) opt->split_factor = 10.;
			if (!opt0.min_chain_weight) opt->min_chain_weight = 40;
			if (strcmp(mode, "pbread1") == 0 || strcmp(mode, "pbread") == 0) {
				opt->flag |= MEM_F_ALL | MEM_F_SELF_OVLP | MEM_F_ALN_REG;
				if (!opt0.max_occ) opt->max_occ = 1000;
				if (!opt0.min_seed_len) opt->min_seed_len = 13;
				if (!opt0.max_chain_extend) opt->max_chain_extend = 25;
				if (opt0.drop_ratio == 0.) opt->drop_ratio = .001;
			} else {
				if (!opt0.min_seed_len) opt->min_seed_len = 17;
				if (!opt0.pen_clip5) opt->pen_clip5 = 0;
				if (!opt0.pen_clip3) opt->pen_clip3 = 0;
			}
		} else {
			fprintf(stderr, "[E::%s] unknown read type '%s'\n", __func__, mode);
			return 1; // FIXME memory leak
		}
	} else update_a(opt, &opt0);
//	if (opt->T < opt->min_HSP_score) opt->T = opt->min_HSP_score; // TODO: tie ->T to MEM_HSP_COEF
	bwa_fill_scmat(opt->a, opt->b, opt->mat);

	gettimeofday (&sTime, NULL);
	if ((idx = bwa_idx_load(argv[optind], BWA_IDX_ALL)) == 0) return 1; // FIXME: memory leak
	gettimeofday (&eTime, NULL);
	idxE=(eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);

	char readIndexFile[500];
	sprintf (readIndexFile, "%s.idx", argv[optind + 1]);
	int status = stat (readIndexFile, &buf);
	assert (status == 0);

        int indFileLen = buf.st_size;
        int numOffsets = indFileLen/8;
	int curInd = rank;
        long *readOffsets = (long *)malloc(numOffsets * sizeof(long));
        FILE *fp_rOfs = fopen(readIndexFile, "r");
	assert(fp_rOfs != NULL);
        fread(readOffsets, sizeof(long), numOffsets, fp_rOfs);
        fclose(fp_rOfs);
	
	status = stat (argv[optind + 1], &buf);
	assert (status == 0);
	long readFileLen=buf.st_size;

	ko = kopen(argv[optind + 1], &fd);
	if (ko == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 1]);
		return 1;
	}
	fp = gzdopen(fd, "r");
	long startOfs = readOffsets[curInd];
	gzseek (fp, (z_off_t)startOfs, SEEK_SET);
        unsigned long curOfs = startOfs;
	ks = kseq_init(fp);
	if (optind + 2 < argc) {
		if (opt->flag&MEM_F_PE) {
			if (bwa_verbose >= 2)
				fprintf(stderr, "[W::%s] when '-p' is in use, the second query file will be ignored.\n", __func__);
		} else {
			ko2 = kopen(argv[optind + 2], &fd2);
			if (ko2 == 0) {
				if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 2]);
				return 1;
			}
			fp2 = gzdopen(fd2, "r");
			gzseek(fp2, curOfs, SEEK_SET);
			ks2 = kseq_init(fp2);
			opt->flag |= MEM_F_PE;
		}
	}
        bntseq_t *bns1 = idx->bns;
        for (i=0; i<bns1->n_seqs; i++)
        {
		bns1->anns[i].numThr=opt->n_threads;
                bns1->anns[i].curFullCount=(long *)malloc(opt->n_threads*sizeof(long));
                bns1->anns[i].curOtherCount=(long *)malloc(opt->n_threads*sizeof(long));
                bns1->anns[i].maxFullCount=(long *)malloc(opt->n_threads*sizeof(long));
                bns1->anns[i].maxOtherCount=(long *)malloc(opt->n_threads*sizeof(long));
		bns1->anns[i].Fr = (fullRec **)malloc(opt->n_threads*sizeof(fullRec *));
		bns1->anns[i].Or = (otherRec **)malloc(opt->n_threads*sizeof(otherRec *));
		for (j=0; j<opt->n_threads; j++)
		{
			bns1->anns[i].curFullCount[j]=0;
			bns1->anns[i].curOtherCount[j]=0;
			bns1->anns[i].maxFullCount[j]=1000;
			bns1->anns[i].maxOtherCount[j]=1000;
	                bns1->anns[i].Fr[j] = (fullRec *)malloc(bns1->anns[i].maxFullCount[j]*sizeof(fullRec));
        	        bns1->anns[i].Or[j] = (otherRec *)malloc(bns1->anns[i].maxOtherCount[j]*sizeof(otherRec));
                	assert (bns1->anns[i].Fr != NULL);
	                assert (bns1->anns[i].Or != NULL);
		}
        }
	FILE *outsam = NULL;

	char hostname[200];
        gethostname(hostname, 200);

	int bytesPerRead=9999;
	struct timeval stTime, etTime;
	gettimeofday (&sTime, NULL);
	while (curInd < numOffsets)
       {
		long ofs1=gztell(fp),ofs2;
		long remainingBytes=(curInd<(numOffsets-1))?(readOffsets[curInd+1]-readOffsets[curInd]):(readFileLen-readOffsets[curInd]);
		long seqsToRead=opt->chunk_size/* * opt->n_threads*/; 
		long totalLen=ofs1;
		while ((seqs = bseq_read(seqsToRead, &n, ks, ks2)) != 0) {
			int64_t size = 0;
			ofs1=gztell(fp);
			if (fp2!=NULL)
	                	ofs2=gztell(fp2);
			if ((opt->flag & MEM_F_PE) && (n&1) == 1) {
				if (bwa_verbose >= 2)
					fprintf(stderr, "[W::%s] odd number of reads in the PE mode; last read dropped\n", __func__);
				n = n>>1<<1;
			}
			for (i = 0; i < n; ++i) size += seqs[i].l_seq;
			if (bwa_verbose >= 3)
				fprintf(stderr, "[M::%s] read %d sequences (%ld bp)...\n", __func__, n, (long)size);

			gettimeofday (&stTime, NULL);
			mem_process_seqs(opt, idx->bwt, idx->bns, idx->pac, n_processed, n, seqs, pes0);
			gettimeofday (&etTime, NULL);
			long elapsed = (etTime.tv_sec * 1000000 + etTime.tv_usec) - (stTime.tv_sec * 1000000 + stTime.tv_usec);
			compE+=elapsed;
			n_processed += n;
			for (i = 0; i < n; ++i) {
				if ((i%2) == 1)
                                {
                                        int nameLen = strlen(seqs[i].name);
					if (seqs[i].comment) nameLen+= strlen(seqs[i].comment);
                                        totalLen+=203+4+2+nameLen;
                                        remainingBytes-=(203+4+2+nameLen);
                                        if (bytesPerRead > (203+4+2+nameLen))
                                                bytesPerRead = 203+4+2+nameLen;
                                        long remSeqs = remainingBytes/bytesPerRead;
                                        if ((remainingBytes%bytesPerRead) > 0) remSeqs++;
                                        remSeqs*=2;
                                        if ((remSeqs*101) < seqsToRead) seqsToRead = remSeqs*101;
                                }
                                if (remainingBytes <= 0)
                                {
                                        break;
                                }
				free(seqs[i].name); 
				free(seqs[i].comment); 
				free(seqs[i].seq); 
				free(seqs[i].qual); 
			}
			ofs1=gztell(fp);
			if (fp2 != NULL)
	                	ofs2=gztell(fp2);
			if (remainingBytes <= 0)
			{
				break;
			}
			free(seqs);
			seqs=NULL;
		}
		curInd = curInd + numtasks;
                startOfs = readOffsets[curInd];
                gzseek (fp, (z_off_t)startOfs, SEEK_SET);
                curOfs = startOfs;
                kseq_destroy(ks);
                ks = kseq_init(fp);
                if (fp2 != NULL)
                {
                        gzseek (fp2, (z_off_t)startOfs, SEEK_SET);
                        kseq_destroy(ks2);
                        ks2 = kseq_init(fp2);
                }
	}
	if (idx->bwt) {bwt_destroy(idx->bwt); idx->bwt=NULL;}
        if (idx->pac) {free(idx->pac); idx->pac=NULL;}

	fprintf (stderr,"Rank %d, align Done\n", rank);
	long *totalRecords = (long *)malloc(bns1->n_seqs*numtasks*2/*aeb,aib*/*sizeof(long));
	gettimeofday (&eTime, NULL);
	long elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
	/* layout first 2*bns1->n_seqs - aeb count for all contigs for PE0 followed by aib of PE0 ....*/
	assert(totalRecords!=NULL);
	bzero (totalRecords, bns1->n_seqs*numtasks*2/*aeb,aib*/*sizeof(long));
	for (i=0; i<bns1->n_seqs; i++)
        {
		totalRecords[(rank*2*bns1->n_seqs)+i]=bns1->anns[i].curFullCount[0];
		totalRecords[(rank*2*bns1->n_seqs)+bns1->n_seqs+i]=bns1->anns[i].curOtherCount[0];
		for (j=1; j<opt->n_threads; j++)
		{
			totalRecords[(rank*2*bns1->n_seqs)+i] += bns1->anns[i].curFullCount[j];
			totalRecords[(rank*2*bns1->n_seqs)+bns1->n_seqs+i] += bns1->anns[i].curOtherCount[j];
		}
		if (opt->n_threads > 1)
		{
			if (totalRecords[(rank*2*bns1->n_seqs)+i] > 0)
			{
				bns1->anns[i].Fr[0]=(fullRec *)realloc(bns1->anns[i].Fr[0], totalRecords[(rank*2*bns1->n_seqs)+i]*sizeof(fullRec));
				assert (bns1->anns[i].Fr[0] != NULL);
			}
			if (totalRecords[(rank*2*bns1->n_seqs)+bns1->n_seqs+i] > 0)
			{
				bns1->anns[i].Or[0]=(otherRec *)realloc(bns1->anns[i].Or[0], totalRecords[(rank*2*bns1->n_seqs)+bns1->n_seqs+i] * sizeof(otherRec));
				assert (bns1->anns[i].Or[0] != NULL);
			}
		}
		for (j=1; j<opt->n_threads; j++)
                {
			memcpy (&(bns1->anns[i].Fr[0][bns1->anns[i].curFullCount[0]]), bns1->anns[i].Fr[j], bns1->anns[i].curFullCount[j]*sizeof(fullRec));
			memcpy (&(bns1->anns[i].Or[0][bns1->anns[i].curOtherCount[0]]), bns1->anns[i].Or[j], bns1->anns[i].curOtherCount[j]*sizeof(otherRec));
			bns1->anns[i].curFullCount[0]+=bns1->anns[i].curFullCount[j];
			bns1->anns[i].curOtherCount[0]+= bns1->anns[i].curOtherCount[j];
			free (bns1->anns[i].Fr[j]);
			bns1->anns[i].Fr[j] = NULL;
			free (bns1->anns[i].Or[j]);
			bns1->anns[i].Or[j] = NULL;
		}
		assert (bns1->anns[i].curFullCount[0] == totalRecords[(rank*2*bns1->n_seqs)+i]);
		assert (bns1->anns[i].curOtherCount[0] == totalRecords[(rank*2*bns1->n_seqs)+bns1->n_seqs+i]);
        }
#ifdef USE_MPI
	long *sendBuf = (long *) malloc (2*bns1->n_seqs*sizeof(long));
	memcpy(sendBuf, &totalRecords[rank*2*bns1->n_seqs], 2*bns1->n_seqs*sizeof(long));
	MPI_Allgather(sendBuf, 2*bns1->n_seqs, MPI_LONG, totalRecords, 2*bns1->n_seqs, MPI_LONG, MPI_COMM_WORLD);
	free(sendBuf);
#endif
	long *reducedTotal = (long *)malloc(bns1->n_seqs * sizeof(long));
	assert (reducedTotal!=NULL);
	bzero (reducedTotal, bns1->n_seqs * sizeof(long));
	long grandTotal=0;
	int index=0;
	for (i=0; i<numtasks; i++)
	{
		for (j=0; j<bns1->n_seqs; j++)
		{
			reducedTotal[j]+=totalRecords[index]+totalRecords[index+bns1->n_seqs];
			grandTotal+=totalRecords[index]+totalRecords[index+bns1->n_seqs];
			index++;
		}
		index+=bns1->n_seqs;
	}
	long workPerProcess = grandTotal/numtasks;
	for (j=1; j<bns1->n_seqs; j++)
		reducedTotal[j]+=reducedTotal[j-1];
	int *owner=(int *)malloc(bns1->n_seqs * sizeof (int));
	assert(owner != NULL);
	int startContig=0, endContig=0;
	owner[0]=0;
	for (i=0; i<numtasks; i++)
	{
		if (i>0)
		{
			startContig=endContig+1;
			endContig=startContig;
			owner[startContig]=i;
		}
		while((endContig < (bns1->n_seqs-1)) && (reducedTotal[endContig+1] < ((i+1)*workPerProcess)))
		{
			endContig++;
			owner[endContig]=i;
		}
	}
	for (i=endContig; i<bns1->n_seqs; i++) owner[i]=numtasks-1;
	for (i=0; i<bns1->n_seqs; i++)
	{
		if (owner[i]==rank)
		{
			startContig=i;
			break;
		}
	}
	for (i=bns1->n_seqs-1; i>=0; i--)
	{
		if (owner[i]==rank)
		{
			endContig=i;
			break;
		}
	}

#ifdef USE_MPI
	MPI_Datatype aebtype,aibtype;
	createAebType (&aebtype);
	createAibType (&aibtype);
	MPI_Type_commit(&aebtype);
	MPI_Type_commit(&aibtype);
#endif

	/* Create and load the arguments to alltoallv */
        fullRec *sendAebBuffer = NULL;
        int sendAebCount = 0;
        otherRec *sendAibBuffer = NULL;
        int sendAibCount = 0;
        fullRec *recvAebBuffer = NULL;
        int recvAebCount = 0;
        otherRec *recvAibBuffer = NULL;
        int recvAibCount = 0;
        int *sendAebCounts = (int *)malloc( numtasks * sizeof(int) );
        int *recvAebCounts = (int *)malloc( numtasks * sizeof(int) );
        int *sendAibCounts = (int *)malloc( numtasks * sizeof(int) );
        int *recvAibCounts = (int *)malloc( numtasks * sizeof(int) );
        int *rAebdispls = (int *)malloc( numtasks * sizeof(int) );
        int *sAebdispls = (int *)malloc( numtasks * sizeof(int) );
        int *rAibdispls = (int *)malloc( numtasks * sizeof(int) );
        int *sAibdispls = (int *)malloc( numtasks * sizeof(int) );
#ifdef USE_MPI
        if (!sendAebCounts || !recvAebCounts || !sendAibCounts || !recvAibCounts || !rAebdispls || !sAebdispls || !rAibdispls || !sAibdispls) {
            fprintf( stderr, "Could not allocate arg items!\n" );fflush(stderr);
            MPI_Abort( MPI_COMM_WORLD, 1 );
        }
#endif
        bzero (sendAebCounts, numtasks * sizeof(int));
        bzero (sendAibCounts, numtasks * sizeof(int));
        bzero (recvAebCounts, numtasks * sizeof(int));
        bzero (recvAibCounts, numtasks * sizeof(int));
        for (j=0; j<numtasks; j++)
        {
                sAebdispls[j]=-1;
                rAebdispls[j]=-1;
                sAibdispls[j]=-1;
                rAibdispls[j]=-1;
        }
        for (i=0; i<bns1->n_seqs; i++)
        {
                if ((i>=startContig) && (i<=endContig))
                {
                        for (j=0; j<numtasks; j++)
                        {
                                if (j!=rank)
                                {
                                        recvAebCount+=totalRecords[(2*bns1->n_seqs*j)+i];
                                        recvAibCount+=totalRecords[(2*bns1->n_seqs*j)+bns1->n_seqs+i];
                                        if (rAebdispls[j]==-1)
                                                rAebdispls[j]=recvAebCounts[j];
                                        if (rAibdispls[j]==-1)
                                                rAibdispls[j]=recvAibCounts[j];
                                        recvAebCounts[j]+=totalRecords[(2*bns1->n_seqs*j)+i];
                                        recvAibCounts[j]+=totalRecords[(2*bns1->n_seqs*j)+bns1->n_seqs+i];
                                }
                        }
                }
                else
                {
                        sendAebCount+=totalRecords[(2*bns1->n_seqs*rank)+i];
                        sendAibCount+=totalRecords[(2*bns1->n_seqs*rank)+bns1->n_seqs+i];
                        if (sAebdispls[owner[i]]==-1)
                                sAebdispls[owner[i]]=sendAebCounts[owner[i]];
                        if (sAibdispls[owner[i]]==-1)
                                sAibdispls[owner[i]]=sendAibCounts[owner[i]];
                        sendAebCounts[owner[i]]+=totalRecords[(2*bns1->n_seqs*rank)+i];
                        sendAibCounts[owner[i]]+=totalRecords[(2*bns1->n_seqs*rank)+bns1->n_seqs+i];
                }
        }
	if (rank > 0)
        {
                rAebdispls[rank]=rAebdispls[rank-1]+recvAebCounts[rank-1];
                rAibdispls[rank]=rAibdispls[rank-1]+recvAibCounts[rank-1];
                sAebdispls[rank]=sAebdispls[rank-1]+sendAebCounts[rank-1];
                sAibdispls[rank]=sAibdispls[rank-1]+sendAibCounts[rank-1];
        }
        else if (rank < (numtasks-1))
        {
                rAebdispls[rank]=rAebdispls[rank+1];
                rAibdispls[rank]=rAibdispls[rank+1];
                sAebdispls[rank]=sAebdispls[rank+1];
                sAibdispls[rank]=sAibdispls[rank+1];
        }
        sendAebBuffer = (fullRec *)malloc(sendAebCount*sizeof(fullRec));
        sendAibBuffer = (otherRec *)malloc(sendAibCount*sizeof(otherRec));
        recvAebBuffer = (fullRec *)malloc(recvAebCount*sizeof(fullRec));
        recvAibBuffer = (otherRec *)malloc(recvAibCount*sizeof(otherRec));
        assert (sendAebBuffer !=NULL);
        assert (sendAibBuffer !=NULL);
        assert (recvAebBuffer !=NULL);
        assert (recvAibBuffer !=NULL);

        long totalAeb=0;
        long totalAib=0;
	for (i=0; i<bns1->n_seqs; i++)
        {
                if ((i>=startContig) && (i<=endContig))
                {
                }
                else
                {
                        if (totalRecords[(2*bns1->n_seqs*rank)+i] > 0)
                        {
                                memcpy (&(sendAebBuffer[totalAeb]), bns1->anns[i].Fr[0], totalRecords[(2*bns1->n_seqs*rank)+i]*sizeof(fullRec));
                                totalAeb+=totalRecords[(2*bns1->n_seqs*rank)+i];
				free(bns1->anns[i].Fr[0]);
				bns1->anns[i].Fr[0]=NULL;
                        }
                        if (totalRecords[(2*bns1->n_seqs*rank)+bns1->n_seqs+i] > 0)
                        {
                                memcpy (&(sendAibBuffer[totalAib]), bns1->anns[i].Or[0], totalRecords[(2*bns1->n_seqs*rank)+bns1->n_seqs+i] * sizeof(otherRec));
                                totalAib+=totalRecords[(2*bns1->n_seqs*rank)+bns1->n_seqs+i];
				free(bns1->anns[i].Or[0]);
				bns1->anns[i].Or[0]=NULL;
                        }
                }
        }
        assert (totalAeb == sendAebCount);
        assert (totalAib == sendAibCount);

#ifdef USE_MPI
        MPI_Alltoallv(sendAebBuffer, sendAebCounts, sAebdispls, aebtype, recvAebBuffer, recvAebCounts, rAebdispls, aebtype, MPI_COMM_WORLD);
        MPI_Alltoallv(sendAibBuffer, sendAibCounts, sAibdispls, aibtype, recvAibBuffer, recvAibCounts, rAibdispls, aibtype, MPI_COMM_WORLD);
#endif

        free(sendAebBuffer);
        free(sendAibBuffer);

        fullRec **allFr = (fullRec **)malloc(bns1->n_seqs*sizeof(fullRec *));
        otherRec **allOr = (otherRec **)malloc (bns1->n_seqs*sizeof(otherRec *));
        long *curOffsFr = (long *)malloc (bns1->n_seqs*sizeof(long));
        long *curOffsOr = (long *)malloc (bns1->n_seqs*sizeof(long));
        assert(curOffsFr!=NULL);
        assert(curOffsOr!=NULL);
        bzero (allFr, bns1->n_seqs*sizeof(fullRec *));
        bzero (allOr, bns1->n_seqs*sizeof(otherRec *));
        bzero (curOffsFr, bns1->n_seqs*sizeof(long));
        bzero (curOffsOr, bns1->n_seqs*sizeof(long));
        for (i=startContig; i<=endContig; i++)
        {
                long totalAeb=0;
                for (j=0; j<numtasks; j++)
                {
                        totalAeb+=totalRecords[(2*bns1->n_seqs*j)+i];
                }
                if (totalAeb > 0)
                {
                        allFr[i] = (fullRec *)malloc(totalAeb * sizeof(fullRec));
                        assert (allFr[i]!=NULL);
                        if (bns1->anns[i].curFullCount[0] > 0)
                        {
                                memcpy(&(allFr[i][curOffsFr[i]]), bns1->anns[i].Fr[0], bns1->anns[i].curFullCount[0]*sizeof(fullRec));
                                curOffsFr[i]+=bns1->anns[i].curFullCount[0];
                        }
                        int curRcvOfs = 0;
                        for (j=0; j<numtasks; j++)
                        {
                                if (j!=rank)
                                {
                                        int k=0;
                                        int contigOffset=curRcvOfs;
                                        for (k=startContig; k<i; k++)
                                        {
                                                contigOffset+=totalRecords[(2*bns1->n_seqs*j)+k];
                                        }
                                        if (totalRecords[(2*bns1->n_seqs*j)+i] > 0)
                                        {
                                                memcpy(&(allFr[i][curOffsFr[i]]), &(recvAebBuffer[contigOffset]), totalRecords[(2*bns1->n_seqs*j)+i]*sizeof(fullRec));
                                                curOffsFr[i]+=totalRecords[(2*bns1->n_seqs*j)+i];
                                        }
                                }
                                curRcvOfs+=recvAebCounts[j];
                        }

                }
                assert (curOffsFr[i]==totalAeb);
                long totalAib=0;
                for (j=0; j<numtasks; j++)
                {
                        totalAib+=totalRecords[(2*bns1->n_seqs*j)+bns1->n_seqs+i];
                }
                if (totalAib > 0)
                {
                        allOr[i] = (otherRec *)malloc(totalAib * sizeof(otherRec));
                        assert (allOr[i]!=NULL);
                        if (bns1->anns[i].curOtherCount[0] > 0)
                        {
                                memcpy(&(allOr[i][curOffsOr[i]]), bns1->anns[i].Or[0], bns1->anns[i].curOtherCount[0]*sizeof(otherRec));
                                curOffsOr[i]+=bns1->anns[i].curOtherCount[0];
                        }
                        int curRcvOfs = 0;
                        for (j=0; j<numtasks; j++)
                        {
                                if (j!=rank)
                                {
                                        int k=0;
                                        int contigOffset=curRcvOfs;
                                        for (k=startContig; k<i; k++)
                                        {
                                                contigOffset+=totalRecords[(2*bns1->n_seqs*j)+bns1->n_seqs+k];
                                        }
                                        if (totalRecords[(2*bns1->n_seqs*j)+bns1->n_seqs+i] > 0)
                                        {
                                                memcpy(&(allOr[i][curOffsOr[i]]), &(recvAibBuffer[contigOffset]), totalRecords[(2*bns1->n_seqs*j)+bns1->n_seqs+i]*sizeof(otherRec));
                                                curOffsOr[i]+=totalRecords[(2*bns1->n_seqs*j)+bns1->n_seqs+i];
                                        }
                                }
                                curRcvOfs+=recvAibCounts[j];
                        }

                }
                assert (curOffsOr[i]==totalAib);
                qsort(allFr[i], curOffsFr[i], sizeof(fullRec), cmpFullRec);
                qsort(allOr[i], curOffsOr[i], sizeof(otherRec), cmpOtherRec);
		int numThr=opt->n_threads;
		parsnip_init(rank, bns1->anns[i].len, argv[optind], out_prefix, bns1->n_seqs, i, numThr);
                parsnip_data *pData=(parsnip_data *)malloc(opt->n_threads*sizeof(parsnip_data));
                pthread_t *thr= (pthread_t *)malloc(opt->n_threads*sizeof(pthread_t));
                assert(pData!=NULL);
                long *nextOfs=(long *)malloc(opt->n_threads*sizeof(long));
                long *remLen=(long *)malloc(opt->n_threads*sizeof(long));
                long *aebStart=(long *)malloc(opt->n_threads*sizeof(long));
                assert((nextOfs != NULL) && (remLen != NULL) && (aebStart !=NULL));
                for (j=0; j<numThr; j++)
                {
                        pthread_attr_t attr;
                        pthread_attr_init(&attr);
                        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
                        aebStart[j]=curOffsFr[i]/numThr*j;
                        long len = curOffsFr[i]/numThr*(j+1);

                        len = (j==(numThr-1))?curOffsFr[i]:len;
                        len -= aebStart[j];
                        nextOfs[j]=len;
                        remLen[j]=0;
                        pData[j].fR=&(allFr[i][aebStart[j]]);
                        pData[j].len=len;
                        pData[j].contig=i;
                        pData[j].tid=j;
                        pData[j].nextOfs = &(nextOfs[j]);
                        pData[j].remLen = &(remLen[j]);
                        pthread_create(&thr[j], &attr, parsnip_aeb, &pData[j]);
                }
		 for (j=0; j<numThr; j++)
                {
                        pthread_join(thr[j], 0);
                        pData[j].fR=&(allFr[i][aebStart[j]+nextOfs[j]]);
                        pData[j].len=remLen[j];
                        pData[j].nextOfs=NULL;
                        pData[j].remLen=NULL;
                        if (pData[j].len > 0)
                                parsnip_aeb(&pData[j]);
                }
		if (curOffsOr[i] > 0)
	                parsnip_aib(allOr[i], curOffsOr[i], i);
		for (j=0; j<opt->n_threads; j++)
                {
                        pthread_attr_t attr;
                        pthread_attr_init(&attr);
                        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
                        long seqStart=bns1->anns[i].len/opt->n_threads*j;
                        long len = bns1->anns[i].len/opt->n_threads*(j+1);

                        len = (j==(opt->n_threads-1))?bns1->anns[i].len:len;
                        len -= seqStart;
                        pData[j].fR=(fullRec *)bns1->anns[i].name;
                        pData[j].len=len;
                        pData[j].contig=(int)seqStart;
                        pData[j].tid=j;
                        pthread_create(&thr[j], &attr, parsnip_write, &pData[j]);
                }
                for (j=0; j<opt->n_threads; j++)
                {
                        pthread_join(thr[j], 0);
                }
                free (pData);
                free (thr);
                free (nextOfs);
                free (remLen);
                free (aebStart);

                gettimeofday (&eTime, NULL);
                elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);

                if (allFr[i] != NULL)
                        free(allFr[i]);
                if (allOr[i] != NULL)
                        free(allOr[i]);

        }
#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
        char suffix[500];
        char vcfBuf[10000000];
        if (rank==0)
        {
                sprintf (suffix, "%s.vcf", out_prefix);
                FILE *fp_vcf = fopen(suffix, "w");
                assert (fp_vcf != NULL);
                fprintf (fp_vcf, "##fileformat=VCFv4.1\n##source=SPRITE3.0\n##reference=%s\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth at the locus\">\n##INFO=<ID=ACTG,Number=4,Type=Integer,Description=\"Count for each allele type\">\n##INFO=<ID=FC,Number=1,Type=Integer,Description=\"Alternate allele count in forward strands\">\n##INFO=<ID=RC,Number=1,Type=Integer,Description=\"Alternate allele count in reverse strands\">\n##INFO=<ID=FT,Number=1,Type=Integer,Description=\"Total forward strands\">\n##INFO=<ID=RT,Number=1,Type=Integer,Description=\"Total reverse strands\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tvcf\n",argv[optind]);
                int nread=0;
                for (i=0; i<bns1->n_seqs; i++)
                {
                        for (j=0; j<opt->n_threads; j++)
                        {
                                sprintf (suffix, "%s_c%d_t%d.vcf", out_prefix, i, j);
                                FILE *tmpVcfFp = fopen (suffix, "r");
                                if (tmpVcfFp != NULL)
                                {
                                        nread = fread(vcfBuf, 1, 10000000, tmpVcfFp);
                                        while (nread > 0)
                                        {
                                                fwrite (vcfBuf, 1, nread, fp_vcf);
                                                nread = fread(vcfBuf, 1, 10000000, tmpVcfFp);
                                        }
                                        fclose (tmpVcfFp);
                                        unlink(suffix);
                                }
                        }
                }
                fclose (fp_vcf);
        }
	free(opt);
	free (allFr);
	free (allOr);
	free(curOffsFr);
	free(curOffsOr);

        free(recvAebBuffer);
        free(recvAibBuffer);

        free(sendAebCounts);
        free(recvAebCounts);
        free(sendAibCounts);
        free(recvAibCounts);
        free(rAebdispls);
        free(sAebdispls);
        free(rAibdispls);
        free(sAibdispls);

	free (reducedTotal);
	for (i=0; i<bns1->n_seqs; i++)
        {
		if (bns1->anns[i].Fr[0] != NULL)
			free(bns1->anns[i].Fr[0]);
                free(bns1->anns[i].Fr);
		if (bns1->anns[i].Or[0] != NULL)
			free(bns1->anns[i].Or[0]);
                free(bns1->anns[i].Or);
		free(bns1->anns[i].curFullCount);
		free(bns1->anns[i].curOtherCount);
		free(bns1->anns[i].maxFullCount);
		free(bns1->anns[i].maxOtherCount);
	}
	free(totalRecords);
	kseq_destroy(ks);
	bwa_idx_destroy(idx);
	err_gzclose(fp); kclose(ko);
	if (ks2) {
		kseq_destroy(ks2);
		err_gzclose(fp2); kclose(ko2);
	}
	if (outsam != NULL)
		fclose (outsam);

	gettimeofday (&eTime, NULL);
	elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
	printf ("Rank %d: Total time %lf s\n", rank, (double)elapsed/1000000);
#ifdef USE_MPI
	MPI_Finalize();
#endif
	return 0;
}
