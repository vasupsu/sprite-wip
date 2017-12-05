#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <assert.h>
#include "minimap.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

#define MM_VERSION "0.2-r124-dirty"

char **refFasta = NULL;

void liftrlimit()
{
#ifdef __linux__
	struct rlimit r;
	getrlimit(RLIMIT_AS, &r);
	r.rlim_cur = r.rlim_max;
	setrlimit(RLIMIT_AS, &r);
#endif
}

void split_fasta (int numChunks, int32_t *len, int numContigs, int *numOutChunks, int *numMaxChunks, int *maxChunkSize)
{
	int32_t maxContigSize = 0;
	int i=0;
	for (i=0; i<numContigs; i++)
	{
		if (len[i] > maxContigSize) maxContigSize = len[i];
	}
	int curChunks = 0, curContigSize=maxContigSize;
	while (curChunks < numChunks)
	{
		curChunks = 0;
		for (i=0; i<numContigs; i++)
		{
			curChunks += (len[i]/curContigSize)+(len[i]%curContigSize > 0);
		}
		printf ("curChunkSize %d --> numChunks %d\n", curContigSize, curChunks);
		if (curChunks < numChunks)
			curContigSize = curContigSize / 2;
	}
	if (curContigSize < maxContigSize)
		curContigSize = curContigSize * 2;
	int maxChunksPerContig = 0;
	curChunks = 0;
	for (i=0; i<numContigs; i++)
	{
		int curChunksPerContig = (len[i]/curContigSize)+(len[i]%curContigSize > 0);
		curChunks += curChunksPerContig;
		if (curChunksPerContig > maxChunksPerContig) maxChunksPerContig = curChunksPerContig;
	}
	*numOutChunks = curChunks;
	*numMaxChunks = maxChunksPerContig;
	*maxChunkSize = curContigSize;
}
void readFasta (char *fileName, int numContigs, int32_t *len)
{
	int i,j;
	struct stat buf;
	int res = stat (fileName, &buf);
	assert (res==0);
	size_t fSize = buf.st_size;
	char *fastaContents = (char *)malloc (fSize);
	assert (fastaContents != NULL);
	FILE *fp = fopen(fileName, "r");
	fread (fastaContents, 1, fSize, fp);
	fclose (fp);
//	printf ("fastaContents[%lu]=%s\n", fSize-100, &fastaContents[fSize-100]);

	size_t curInd=0;
	for (i=0; i<numContigs; i++)
	{
		printf ("%d - ", i);
		while ((curInd < fSize) && (fastaContents[curInd]!='\n')) 
		{
			printf ("%c", fastaContents[curInd]);
			curInd++;
		}
		printf ("\n");
		assert (curInd < fSize);
		curInd++;
		for (j=0; j<len[i]; j++)
		{
			refFasta[i][j]=fastaContents[curInd++];
/*			if ((contigs[i].len-j)<100)
				printf ("%c", refFasta[i][j]);*/
			if (refFasta[i][j]=='\n')
				j--;
			else if (refFasta[i][j] > 96)
				refFasta[i][j] = refFasta[i][j] - 32;
		}
//		printf ("%s\n", &(refFasta[i][contigs[i].len-100]));
		assert (fastaContents[curInd] == '\n');
		curInd++;
	}
	printf ("curInd %lu fSize %lu\n", curInd, fSize);
	free (fastaContents);
}
int main(int argc, char *argv[])
{
	mm_mapopt_t opt;
	int i, c, k = 15, w = -1, b = MM_IDX_DEF_B, n_threads = 3, keep_name = 1, is_idx = 0;
	int tbatch_size = 100000000;
	int numTasks = 16, rank = 0;
#ifdef USE_MPI
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numTasks);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	char hostname[100];
	gethostname (hostname, 100);
	printf ("numTasks %d rank %d hostname %s\n", numTasks, rank, hostname);
#endif
	uint64_t ibatch_size = 4000000000ULL;
	float f = 0.001;
	bseq_file_t *fp = 0;
	char *fnw = 0, *oprefix=NULL;
	FILE *fpr = 0, *fpw = 0;

	liftrlimit();
	mm_realtime0 = realtime();
	printf ("Start time %lf\n", mm_realtime0);
	mm_mapopt_init(&opt);

	while ((c = getopt(argc, argv, "o:w:k:B:b:t:r:c:f:Vv:NOg:I:d:lRPSpT:m:L:Dx:")) >= 0) {
		if (c == 'w') w = atoi(optarg);
		else if (c == 'k') k = atoi(optarg);
		else if (c == 'b') b = atoi(optarg);
		else if (c == 'r') opt.radius = atoi(optarg);
		else if (c == 'c') opt.min_cnt = atoi(optarg);
		else if (c == 'm') opt.merge_frac = atof(optarg);
		else if (c == 'f') f = atof(optarg);
		else if (c == 't') n_threads = atoi(optarg);
		else if (c == 'v') mm_verbose = atoi(optarg);
		else if (c == 'g') opt.max_gap = atoi(optarg);
		else if (c == 'N') keep_name = 0;
		else if (c == 'd') fnw = optarg;
		else if (c == 'o') oprefix = optarg;//vasu
		else if (c == 'l') is_idx = 1;
		else if (c == 'R') opt.flag |= MM_F_WITH_REP;
		else if (c == 'P') opt.flag &= ~MM_F_WITH_REP;
		else if (c == 'D') opt.flag |= MM_F_NO_SELF;
		else if (c == 'O') opt.flag |= MM_F_NO_ISO;
		else if (c == 'S') opt.flag |= MM_F_AVA | MM_F_NO_SELF;
		else if (c == 'p') opt.flag |= MM_F_PE;//vasu
		else if (c == 'T') opt.sdust_thres = atoi(optarg);
		else if (c == 'L') opt.min_match = atoi(optarg);
		else if (c == 'V') {
			puts(MM_VERSION);
			return 0;
		} else if (c == 'B' || c == 'I') {
			double x;
			char *p;
			x = strtod(optarg, &p);
			if (*p == 'G' || *p == 'g') x *= 1e9;
			else if (*p == 'M' || *p == 'm') x *= 1e6;
			else if (*p == 'K' || *p == 'k') x *= 1e3;
			if (c == 'B') tbatch_size = (uint64_t)(x + .499);
			else ibatch_size = (uint64_t)(x + .499);
		} else if (c == 'x') {
			if (strcmp(optarg, "ava10k") == 0) {
				opt.flag |= MM_F_AVA | MM_F_NO_SELF;
				opt.min_match = 100;
				opt.merge_frac = 0.0;
				w = 5;
			}
		}
	}
	if (w < 0) w = (int)(.6666667 * k + .499);

	if (argc == optind) {
		fprintf(stderr, "Usage: minimap [options] <target.fa> [query.fa] [...]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  Indexing:\n");
		fprintf(stderr, "    -k INT     k-mer size [%d]\n", k);
		fprintf(stderr, "    -w INT     minizer window size [{-k}*2/3]\n");
		fprintf(stderr, "    -I NUM     split index for every ~NUM input bases [4G]\n");
		fprintf(stderr, "    -d FILE    dump index to FILE []\n");
		fprintf(stderr, "    -o PREFIX  output prefix path for AEB and Unmapped output files\n");
		fprintf(stderr, "    -l         index file <target.fa>.mmi is used(overriding -k, -w and -I)\n");
//		fprintf(stderr, "    -b INT     bucket bits [%d]\n", b); // most users would care about this
		fprintf(stderr, "  Mapping:\n");
		fprintf(stderr, "    -f FLOAT   filter out top FLOAT fraction of repetitive minimizers [%.3f]\n", f);
		fprintf(stderr, "    -r INT     bandwidth [%d]\n", opt.radius);
		fprintf(stderr, "    -m FLOAT   merge two chains if FLOAT fraction of minimizers are shared [%.2f]\n", opt.merge_frac);
		fprintf(stderr, "    -c INT     retain a mapping if it consists of >=INT minimizers [%d]\n", opt.min_cnt);
		fprintf(stderr, "    -L INT     min matching length [%d]\n", opt.min_match);
		fprintf(stderr, "    -g INT     split a mapping if there is a gap longer than INT [%d]\n", opt.max_gap);
		fprintf(stderr, "    -T INT     SDUST threshold; 0 to disable SDUST [%d]\n", opt.sdust_thres);
//		fprintf(stderr, "    -D         skip self mappings but keep dual mappings\n"); // too confusing to expose to end users
		fprintf(stderr, "    -S         skip self and dual mappings\n");
		fprintf(stderr, "    -p         paired-end mapping a_1.fq a_2.fq b_1.fq b_2.fq...\n"); //vasu
		fprintf(stderr, "    -O         drop isolated hits before chaining (EXPERIMENTAL)\n");
		fprintf(stderr, "    -P         filtering potential repeats after mapping (EXPERIMENTAL)\n");
//		fprintf(stderr, "    -R         skip post-mapping repeat filtering\n"); // deprecated option for backward compatibility
		fprintf(stderr, "    -x STR     preset (recommended to be applied before other options) []\n");
		fprintf(stderr, "               ava10k: -Sw5 -L100 -m0 (PacBio/ONT all-vs-all read mapping)\n");
		fprintf(stderr, "  Input/Output:\n");
		fprintf(stderr, "    -t INT     number of threads [%d]\n", n_threads);
//		fprintf(stderr, "    -B NUM     process ~NUM bp in each batch [100M]\n");
//		fprintf(stderr, "    -v INT     verbose level [%d]\n", mm_verbose);
//		fprintf(stderr, "    -N         use integer as target names\n");
		fprintf(stderr, "    -V         show version number\n");
		fprintf(stderr, "\nSee minimap.1 for detailed description of the command-line options.\n");
		return 1;
	}
	char idxFileName [100];
	sprintf (idxFileName, "%s.mmi", argv[optind]);
	if (is_idx) fpr = fopen(idxFileName, "rb");
	else fp = bseq_open(argv[optind]);
	if (fnw) fpw = fopen(fnw, "wb");
	double indexTime=0, fileStartTime=0;
	for (;;) {
		mm_idx_t *mi = 0;
		if (fpr) mi = mm_idx_load(fpr);
		else if (!bseq_eof(fp))
			mi = mm_idx_gen(fp, w, k, b, tbatch_size, n_threads, ibatch_size, keep_name);
		if (mi == 0) break;
		if (mm_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built the index for %d target sequence(s)\n",
					__func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), mi->n);
		mm_idx_set_max_occ(mi, f);
		if (mm_verbose >= 3)
			fprintf(stderr, "[M::%s] max occurrences of a minimizer to consider: %d\n", __func__, mi->max_occ);
		if (fpw) mm_idx_dump(fpw, mi);
		refFasta = (char **) malloc(mi->n * sizeof(char *));
		assert (refFasta != NULL);
		for (i=0; i<mi->n; i++)
		{
//			printf ("len[%d]=%d\n", i, mi->len[i]);
			refFasta[i] = (char *)malloc(mi->len[i] + 1);
			assert (refFasta[i] != NULL);
		}
		int numOutChunks=0, numMaxChunks=0, maxChunkSize=0;
		split_fasta (500, mi->len, mi->n, &numOutChunks, &numMaxChunks, &maxChunkSize);
		printf ("MaxChunkSize %d, numOutChunks %d, max num of chunks per contig %d\n", maxChunkSize, numOutChunks, numMaxChunks);
		readFasta (argv[optind], mi->n, mi->len);
		indexTime = realtime();
		printf ("index End Timestamp %lf\n",  indexTime);
		FILE **aebFp = (FILE **)calloc ((numMaxChunks * mi->n + 1), sizeof (FILE *));
		assert (aebFp != NULL);
		if (opt.flag & MM_F_PE)
		{
//			fprintf (stderr, "i %d argc %d\n", optind + 1, argc);
			printf ("main: refFasta %p mi %p\n", refFasta, mi);
			char fileName[500];
			for (i = optind + 1; i < argc; i+=2)//vasu
			{
				sprintf (fileName, "%s.idx4", argv[i]);
				struct stat buf;
				int res = stat (fileName, &buf);
				assert (res==0);
				int numRecs = buf.st_size/sizeof(fileIndex);
				if (numTasks > 1)
					assert ((numRecs > numTasks*n_threads) && ((numRecs % (numTasks*n_threads))==1));
				int chunksPerTask = (numRecs-1)/numTasks;
				int startChunk = rank*chunksPerTask;
				int endChunk = (rank+1)*chunksPerTask;
				fileIndex *fI = (fileIndex *)malloc(numRecs * sizeof(fileIndex));
				assert (fI != NULL);

				FILE *fpFI = fopen (fileName, "r");
				fread (fI, sizeof(fileIndex), numRecs, fpFI);
				fclose (fpFI);
				sprintf (fileName, "%s_%d", oprefix, rank);
				printf ("mm_map_file start time %lf\n", realtime());
				mm_map_file(numTasks, rank, numRecs, numOutChunks, maxChunkSize, numMaxChunks, fI, aebFp, fileName, refFasta, mi, argv[i], argv[i+1], &opt, n_threads, tbatch_size);
				printf ("mm_map_file end time %lf\n", realtime());

				free (fI);
			}
		}
		else
		{
			char fileName[500];
			for (i = optind + 1; i < argc; ++i)
			{
				sprintf (fileName, "%s.idx4", argv[i]);
				struct stat buf;
				int res = stat (fileName, &buf);
				assert (res==0);
				int numRecs = buf.st_size/sizeof(fileIndex);
				if (numTasks > 1)
					assert ((numRecs > numTasks*n_threads) && ((numRecs % (numTasks*n_threads))==1));
				int chunksPerTask = (numRecs-1)/numTasks;
				int startChunk = rank*chunksPerTask;
				int endChunk = (rank+1)*chunksPerTask;
				fileIndex *fI = (fileIndex *)malloc(numRecs * sizeof(fileIndex));
				assert (fI != NULL);
				FILE *fpFI = fopen (fileName, "r");
				fread (fI, sizeof(fileIndex), numRecs, fpFI);
				fclose (fpFI);
				sprintf (fileName, "%s_%d", oprefix, rank);
				mm_map_file(numTasks, rank, numRecs, numOutChunks, maxChunkSize, numMaxChunks, fI, aebFp, fileName, refFasta, mi, argv[i], NULL, &opt, n_threads, tbatch_size);
				free (fI);
			}
		}
		int j=0;
		fileStartTime = realtime();
		printf ("file close start time %lf\n", fileStartTime);
		for (i=0; i<mi->n; i++)
		{
			free (refFasta[i]);
			refFasta[i] = NULL;
			for (j=0; j<numMaxChunks; j++)
				if (aebFp[i*numMaxChunks+j] != NULL) 
				{
//					printf ("aebFp[C%d S%d]=%p\n", i, j, aebFp[i*numMaxChunks+j]);
					
					fclose (aebFp[i*numMaxChunks+j]);
					aebFp[i*numMaxChunks+j] = NULL;
				}
		}
		if (aebFp[mi->n*numMaxChunks] != NULL) 
		{
			fclose (aebFp[mi->n*numMaxChunks]);
			aebFp[mi->n*numMaxChunks] = NULL;
		}
//		printf ("File close time %.3lf sec\n", realtime()-fileStartTime);
		mm_idx_destroy(mi);
		free (aebFp);
		free (refFasta);
	}
	if (fpw) fclose(fpw);
	if (fpr) fclose(fpr);
	if (fp)  bseq_close(fp);

	fprintf(stderr, "[M::%s] Version: %s\n", __func__, MM_VERSION);
	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	printf ("app end time %lf\n", realtime());
	fprintf(stderr, "\n[M::%s] Rank %d Real time: %.3f sec(IndexTime %.3lf sec); FileCloseTime %.3lf sec CPU: %.3f sec\n", __func__, rank, realtime() - mm_realtime0, indexTime-mm_realtime0, realtime()-fileStartTime, cputime());
#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0)
	{
		int pid = getpid();
		char line[2048];
		sprintf (line, "/proc/%d/status", pid);
		FILE *statFile = fopen(line, "r");
		assert (statFile != NULL);
		fgets (line, 2048, statFile);
		while (!feof (statFile))
		{
			if (strstr(line,"VmPeak") || strstr(line,"VmHWM"))
			{
				printf ("[%d] %s", rank, line);
			}
			fgets (line, 2048, statFile);
		}
		fclose (statFile);
	}
	MPI_Finalize();
#endif
	return 0;
}
