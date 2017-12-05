#include <zlib.h>
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
#include "kstring.h"
#include "kvec.h"
#include "utils.h"
#include "kseq.h"
#include "utils.h"
#include "faiTrie.h"
#ifdef USE_MPI
#include <mpi.h>
#endif
KSEQ_DECLARE(gzFile)
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.7.12-r1039"
#endif

extern unsigned char nst_nt4_table[256];

void *kopen(const char *fn, int *_fd);
int kclose(void *a);
char *bwa_pg;
struct stat buf;

typedef struct
{
        uint32_t pos;
        uint8_t seq[30];
        uint16_t flag;
        uint8_t qual;
        uint8_t  matchLen;
}fullRec;

typedef struct
{
        uint32_t pos;
        uint16_t flag;
        uint8_t seq[30];
        uint8_t qual;
        uint8_t n_cigar;
        uint16_t cigar[10];
}otherRec;

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

TrieNode *constructTrie (char *fastaFile)
{
	TrieNode *fai_trie = newTrieNode(' ', -1, -1);
	char faiFile[300];
	char line[200];
	char contigName[100];
        int length=0;
	sprintf (faiFile, "%s.fai", fastaFile);
	FILE *fp = fopen(faiFile, "r");
	assert (fp != NULL);
	fgets(line, 200, fp);
	int contigNo=0;
	while (!feof(fp))
	{
		sscanf(line, "%s\t%d", contigName, &length);
		insertStr (fai_trie, contigName, length, contigNo++);
		fgets(line, 200, fp);
	}
	fclose(fp);
	return fai_trie;
}

uint16_t getNextCigar (char **cigarStr)
{
	int type=-1, ofs=0;
	uint16_t len=0;
	char lenStr[10];
	if (**cigarStr=='\0') return 99;

	while (type==-1)
	{
		switch(**cigarStr)
		{
			case 'M':
				type=0;
				break;
			case 'I':
				type=1;
				break;
			case 'D':
				type=2;
				break;
			case 'S':
				type=3;
				break;
			case 'H':
				type=4;
				break;
			default:
                                if (isdigit(**cigarStr))
                                {
                                        lenStr[ofs++]=**cigarStr;
                                }
				else
				{
					return -1;
				}
		}
		(*cigarStr)++;
	}
	lenStr[ofs]='\0';
	len=(uint16_t)atoi(lenStr);
	return ((len<<4)|type);
}
//return numCigarOps and cigarOpList
int getCigars(char *cigarStr, uint16_t *cigarOpList)
{
	int numOps=0;
	if (*cigarStr=='*') return 0;
	uint16_t cigarOp = getNextCigar (&cigarStr);
	while (cigarOp != 99)
	{
		cigarOpList[numOps++]=cigarOp;
		cigarOp = getNextCigar (&cigarStr);
	}
	return numOps;
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
int main(int argc, char *argv[])
{
	mem_opt_t *opt, opt0;
	int fd, fd2, i, c, n, copy_comment=0, j;
	gzFile fp, fp2 = 0;
	kseq_t *ks, *ks2 = 0;
	bseq1_t *seqs;
	bwaidx_t *idx;
	char *p, *rg_line = 0;
	const char *mode = 0;
	void *ko = 0, *ko2 = 0;
	int numtasks=1, rank=0;
	int64_t n_processed = 0;
	int create_sam = 0;
	mem_pestat_t pes[4], *pes0 = 0;
	kstring_t pg = {0,0,0};
        ksprintf(&pg, "@PG\tID:bwa\tPN:bwa\tVN:%s\tCL:%s", PACKAGE_VERSION, argv[0]);
        for (i = 1; i < argc; ++i) ksprintf(&pg, " %s", argv[i]);
        bwa_pg = pg.s;
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
	while ((c = getopt(argc, argv, "epaFMCSPHYk:c:v:s:r:t:R:A:B:O:E:U:w:L:d:T:Q:D:m:I:N:W:x:G:h:o:")) >= 0) {
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
		else if (c == 's') {create_sam = atoi(optarg); if (create_sam!=1) create_sam=0;}
//		else if (c == 's') opt->split_width = atoi(optarg), opt0.split_width = 1;
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
	if (optind + 1 >= argc || optind + 3 < argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]\n\n");
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
		fprintf(stderr, "       -s INT        create SAM file output (-s 1 to create SAM)\n");
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
	FILE *outsam = NULL;
	if (create_sam)
	{
		char outfile[200];
		if (out_prefix != NULL)
			sprintf (outfile, "%s_%d.sam", out_prefix, rank);
		else
			sprintf (outfile, "%s_%d.sam", argv[optind + 1], rank);
		printf ("outfile %s\n", outfile);
		outsam = fopen(outfile, "w");
		assert (outsam!=NULL);
	        if (!(opt->flag & MEM_F_ALN_REG))
        		bwa_print_sam_hdr2(idx->bns, rg_line, outsam);
	}
	
	char fName[500];
	TrieNode *fai_trie = constructTrie(argv[optind]);

        bntseq_t *bns1 = idx->bns;
	contigInfo *auxInfo = (contigInfo *)malloc(bns1->n_seqs * sizeof(contigInfo));
	assert(auxInfo != NULL);
	fprintf (stderr, "FullRec %lu OtherRec %lu\n", sizeof(fullRec), sizeof(otherRec));
        for (i=0; i<bns1->n_seqs; i++)
        {
		if (out_prefix != NULL)
                	sprintf (fName, "%s_%s_%d.aeb", out_prefix, bns1->anns[i].name, rank);
		else
                	sprintf (fName, "%s_%s_%d.aeb", argv[optind + 1], bns1->anns[i].name, rank);
		auxInfo[i].fullFileName = (char *)malloc(strlen(fName)+1);
		strcpy(auxInfo[i].fullFileName, fName);
                auxInfo[i].fp_fullMatch = NULL/*fopen(fName, "w")*/;
		if (out_prefix != NULL)
                	sprintf (fName, "%s_%s_%d.aib", out_prefix, bns1->anns[i].name, rank);
		else
                	sprintf (fName, "%s_%s_%d.aib", argv[optind + 1], bns1->anns[i].name, rank);
                auxInfo[i].otherFileName = (char *)malloc(strlen(fName)+1);
                strcpy(auxInfo[i].otherFileName, fName);
                auxInfo[i].fp_others = NULL/*fopen(fName, "w")*/;
                auxInfo[i].curFullCount=0;
                auxInfo[i].curOtherCount=0;
                auxInfo[i].outputChunkSize=100000;
                auxInfo[i].Fr = (fullRec *)malloc(auxInfo[i].outputChunkSize*sizeof(fullRec));
                auxInfo[i].Or = (otherRec *)malloc(auxInfo[i].outputChunkSize*sizeof(otherRec));
                assert (auxInfo[i].Fr != NULL);
                assert (auxInfo[i].Or != NULL);
	}

	char hostname[200];
        gethostname(hostname, 200);

	int readLength=0;
	uint32_t alnPos=0;
	int qual=0;
	uint16_t cigarOpList[50];
	uint8_t numCigars=0;
	int contigNo=0;
	uint16_t flag=0;
	fullRec FR;
	otherRec OR;
	int bytesPerRead=9999;
	struct timeval stTime, etTime;
	gettimeofday (&sTime, NULL);
	while (curInd < numOffsets)
        {
		long ofs1=gztell(fp);
		printf ("Process %d Chunk %d/%d\n", rank, curInd, numOffsets);
		long remainingBytes=(curInd<(numOffsets-1))?(readOffsets[curInd+1]-readOffsets[curInd]):(readFileLen-readOffsets[curInd]);
		long seqsToRead=opt->chunk_size/* * opt->n_threads*/; 
		long totalLen=ofs1,ofs2=0;
		while ((seqs = bseq_read(seqsToRead, &n, ks, ks2)) != 0) {
//			fprintf (stderr, "read size %ld, n %ld, bytesPerRead %ld, remBytes %ld\n", seqsToRead, n, bytesPerRead, remainingBytes);
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
				gettimeofday (&stTime, NULL);
				if (seqs[i].sam && create_sam) { err_fputs(seqs[i].sam, outsam); }

				gettimeofday (&etTime, NULL);
				elapsed = (etTime.tv_sec * 1000000 + etTime.tv_usec) - (stTime.tv_sec * 1000000 + stTime.tv_usec);
				writeE += elapsed;

				int token_no=0;
			        char *tok=seqs[i].sam;
			        int len=strlen(tok);
				int seqInd=0;
				FILE *out = NULL;
				int fileType=-1;
				gettimeofday (&stTime, NULL);
				char *curReadId=NULL;
			        for (j=0; j<len; j++)
			        {
			                if ((seqs[i].sam[j]=='\t')||(seqs[i].sam[j]=='\n'))
			                {
			                        seqs[i].sam[j] = '\0';
						if ((curReadId != NULL) && (strcmp(curReadId, tok)==0))
							token_no=0;
						switch (token_no)
						{
							case 0://read name
								curReadId=tok;
								break;
							case 1://flag
								flag = (uint16_t)atoi(tok);
								break;
							case 2://contigName
								contigNo=getContigNo(fai_trie, tok, 0);
								if (contigNo == -1)
									j=len;
								break;
							case 3://pos
								alnPos = (uint32_t)atoi(tok);
								break;
							case 4://qual
								qual = atoi (tok);
								break;
							case 5://cigar
								numCigars = getCigars(tok, cigarOpList);
								if ((numCigars == 1) && ((cigarOpList[0]&0XF)==0))
								{
									fileType=0;//full
									if (auxInfo[contigNo].fp_fullMatch==NULL)
									{
										auxInfo[contigNo].fp_fullMatch = fopen(auxInfo[contigNo].fullFileName, "w");
									}
									out=auxInfo[contigNo].fp_fullMatch;
								}
								else
								{
									if (auxInfo[contigNo].fp_others==NULL)
									{
										auxInfo[contigNo].fp_others = fopen(auxInfo[contigNo].otherFileName, "w");
									}
									out=auxInfo[contigNo].fp_others;
								}
								break;
							case 9://seq
								if (readLength==0)
								{
									readLength=strlen(tok);
								}
								if (fileType==0)
								{
									for (seqInd=0; seqInd<readLength; seqInd++)
									{
										int curInd=seqInd/4;
										char ch=tok[seqInd];
										switch (seqInd&3)
						                                {
						                                        case 0:
						                                           FR.seq[curInd]=(ch<<5)&0xc0;
						                                           break;
						                                        case 1:
						                                           FR.seq[curInd] |= ((ch<<3)&0x30);
                                                					   break;
                                        						case 2:
                                                					   FR.seq[curInd] |= ((ch<<1)&0xc);
                                                					   break;
                                        						case 3:
                                                					   FR.seq[curInd] |= ((ch>>1) & 0x3);
                                                					   break;
                                						}
									}
								}
								else
								{
									for (seqInd=0; seqInd<readLength; seqInd++)
									{
										int curInd=seqInd/4;
										char ch=tok[seqInd];
                                                                                switch (seqInd&3)
                                                                                {
                                                                                        case 0:
                                                                                           OR.seq[curInd]=(ch<<5)&0xc0;
                                                                                           break;
                                                                                        case 1:
                                                                                           OR.seq[curInd] |= ((ch<<3)&0x30);
                                                                                           break;
                                                                                        case 2:
                                                                                           OR.seq[curInd] |= ((ch<<1)&0xc);
                                                                                           break;
                                                                                        case 3:
                                                                                           OR.seq[curInd] |= ((ch>>1) & 0x3);
                                                                                           break;
                                                                                }
									}
								}
								break;
						}
                        			token_no++;
			                        tok=&seqs[i].sam[j+1];
                			}
			        }
				gettimeofday (&etTime, NULL);
				elapsed = (etTime.tv_sec * 1000000 + etTime.tv_usec) - (stTime.tv_sec * 1000000 + stTime.tv_usec);
				overE += elapsed;
				if ((i%2) == 1)
				{
					int nameLen = (seqs[i].comment==NULL)?strlen(seqs[i].name):(strlen(seqs[i].name) + strlen(seqs[i].comment));
					totalLen+=readLength*2+1+4+2+nameLen;
					remainingBytes-=(readLength*2+1+4+2+nameLen);
					if (bytesPerRead > (readLength*2+1+4+2+nameLen))
						bytesPerRead = readLength*2+1+4+2+nameLen;
					long remSeqs = remainingBytes/bytesPerRead;
					if ((remainingBytes%bytesPerRead) > 0) remSeqs++;
					remSeqs*=2;
					if ((remSeqs*readLength) < seqsToRead) seqsToRead = remSeqs*readLength;
				}
//				fprintf (stderr, "remainingBytes %ld\n", remainingBytes);
				if (remainingBytes <= 0)
				{
					break;
				}
				free(seqs[i].name); free(seqs[i].comment); free(seqs[i].seq); free(seqs[i].qual); free(seqs[i].sam);
				if (contigNo == -1) continue;
				if (fileType == 0)
				{
					FR.pos = alnPos;
					FR.flag=flag;
					FR.qual = qual;
					FR.matchLen = readLength;
					gettimeofday (&stTime, NULL);
					if (auxInfo[contigNo].curFullCount == auxInfo[contigNo].outputChunkSize)
					{
						fwrite(auxInfo[contigNo].Fr, sizeof(fullRec), auxInfo[contigNo].curFullCount, auxInfo[contigNo].fp_fullMatch);
						auxInfo[contigNo].curFullCount=0;
					}
					gettimeofday (&etTime, NULL);
					elapsed = (etTime.tv_sec * 1000000 + etTime.tv_usec) - (stTime.tv_sec * 1000000 + stTime.tv_usec);
					writeE += elapsed;
					memcpy(&(auxInfo[contigNo].Fr[auxInfo[contigNo].curFullCount]), &FR, sizeof(fullRec));
					auxInfo[contigNo].curFullCount++;
				}
				else if (numCigars > 0)
				{
					OR.pos = alnPos;
					OR.flag=flag;
					OR.qual = qual;
					OR.n_cigar = (uint8_t)numCigars;
					memcpy (OR.cigar, cigarOpList, 10*sizeof(uint16_t));
					gettimeofday (&stTime, NULL);
					if (auxInfo[contigNo].curOtherCount == auxInfo[contigNo].outputChunkSize)
					{
						fwrite(auxInfo[contigNo].Or, sizeof(otherRec), auxInfo[contigNo].curOtherCount, auxInfo[contigNo].fp_others);
						auxInfo[contigNo].curOtherCount=0;
					}
					gettimeofday (&etTime, NULL);
					elapsed = (etTime.tv_sec * 1000000 + etTime.tv_usec) - (stTime.tv_sec * 1000000 + stTime.tv_usec);
					writeE += elapsed;
					memcpy(&(auxInfo[contigNo].Or[auxInfo[contigNo].curOtherCount]), &OR, sizeof(otherRec));
					auxInfo[contigNo].curOtherCount++;
				}
			}
			ofs1=gztell(fp);
			if (fp2 != NULL)
	                	ofs2=gztell(fp2);
			
			free(seqs);
			seqs=NULL;
			if (remainingBytes <= 0)
			{
				break;
			}
		}
		curInd = curInd + numtasks;
//		printf ("curInd %d\n", curInd);
//		if (curInd > 0) break;
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
#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	free(opt);
	for (i=0; i<bns1->n_seqs; i++)
        {
                free(auxInfo[i].fullFileName);
		if (auxInfo[i].fp_fullMatch != NULL)
		{
			gettimeofday (&stTime, NULL);
			if (auxInfo[i].curFullCount > 0)
			{
				fwrite(auxInfo[i].Fr, sizeof(fullRec), auxInfo[i].curFullCount, auxInfo[i].fp_fullMatch);
			}
			gettimeofday (&etTime, NULL);
			long elapsed = (etTime.tv_sec * 1000000 + etTime.tv_usec) - (stTime.tv_sec * 1000000 + stTime.tv_usec);
			writeE += elapsed;
			fclose (auxInfo[i].fp_fullMatch);
		}
                free(auxInfo[i].otherFileName);
		if (auxInfo[i].fp_others != NULL)
		{
			gettimeofday (&stTime, NULL);
			if (auxInfo[i].curOtherCount > 0)
			{
				fwrite(auxInfo[i].Or, sizeof(otherRec), auxInfo[i].curOtherCount, auxInfo[i].fp_others);
			}
			gettimeofday (&etTime, NULL);
			long elapsed = (etTime.tv_sec * 1000000 + etTime.tv_usec) - (stTime.tv_sec * 1000000 + stTime.tv_usec);
			writeE += elapsed;
                	fclose (auxInfo[i].fp_others);
		}
                free(auxInfo[i].Fr);
                free(auxInfo[i].Or);
        }
	bwa_idx_destroy(idx);
	kseq_destroy(ks);
	err_gzclose(fp); kclose(ko);
	if (ks2) {
		kseq_destroy(ks2);
		err_gzclose(fp2); kclose(ko2);
	}
	free(auxInfo);
	if (outsam != NULL)
		fclose (outsam);
	gettimeofday (&eTime, NULL);
	long elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
	printf ("Rank %d: Comp %lf s, Read %lf s, Write %lf s, bwtRead %lf s, sam2aeb %lf s, Total time %lf s\n", rank, (double)compE/1000000, (double)readE/1000000, (double)writeE/1000000, (double)idxE/1000000, (double)overE/1000000, (double)elapsed/1000000);
#ifdef USE_MPI
	MPI_Finalize();
#endif
	return 0;
}

int main_fastmap(int argc, char *argv[])
{
	int c, i, min_iwidth = 20, min_len = 17, print_seq = 0;
	kseq_t *seq;
	bwtint_t k;
	gzFile fp;
	smem_i *itr;
	const bwtintv_v *a;
	bwaidx_t *idx;

	while ((c = getopt(argc, argv, "w:l:p")) >= 0) {
		switch (c) {
			case 'p': print_seq = 1; break;
			case 'w': min_iwidth = atoi(optarg); break;
			case 'l': min_len = atoi(optarg); break;
		    default: return 1;
		}
	}
	if (optind + 1 >= argc) {
		fprintf(stderr, "Usage: bwa fastmap [-p] [-l minLen=%d] [-w maxSaSize=%d] <idxbase> <in.fq>\n", min_len, min_iwidth);
		return 1;
	}

	fp = xzopen(argv[optind + 1], "r");
	seq = kseq_init(fp);
	if ((idx = bwa_idx_load(argv[optind], BWA_IDX_BWT|BWA_IDX_BNS)) == 0) return 1;
	itr = smem_itr_init(idx->bwt);
	while (kseq_read(seq) >= 0) {
		err_printf("SQ\t%s\t%ld", seq->name.s, seq->seq.l);
		if (print_seq) {
			err_putchar('\t');
			err_puts(seq->seq.s);
		} else err_putchar('\n');
		for (i = 0; i < seq->seq.l; ++i)
			seq->seq.s[i] = nst_nt4_table[(int)seq->seq.s[i]];
		smem_set_query(itr, seq->seq.l, (uint8_t*)seq->seq.s);
		while ((a = smem_next(itr)) != 0) {
			for (i = 0; i < a->n; ++i) {
				bwtintv_t *p = &a->a[i];
				if ((uint32_t)p->info - (p->info>>32) < min_len) continue;
				err_printf("EM\t%d\t%d\t%ld", (uint32_t)(p->info>>32), (uint32_t)p->info, (long)p->x[2]);
				if (p->x[2] <= min_iwidth) {
					for (k = 0; k < p->x[2]; ++k) {
						bwtint_t pos;
						int len, is_rev, ref_id;
						len  = (uint32_t)p->info - (p->info>>32);
						pos = bns_depos(idx->bns, bwt_sa(idx->bwt, p->x[0] + k), &is_rev);
						if (is_rev) pos -= len - 1;
						bns_cnt_ambi(idx->bns, pos, len, &ref_id);
						err_printf("\t%s:%c%ld", idx->bns->anns[ref_id].name, "+-"[is_rev], (long)(pos - idx->bns->anns[ref_id].offset) + 1);
					}
				} else err_puts("\t*");
				err_putchar('\n');
			}
		}
		err_puts("//");
	}

	smem_itr_destroy(itr);
	bwa_idx_destroy(idx);
	kseq_destroy(seq);
	err_gzclose(fp);
	return 0;
}
