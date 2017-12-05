#include "sam.h"
#include "hts.h"
#include "bgzf.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <omp.h>

#define SEGMENT_SIZE 1000000
typedef struct {
	int contigNo;
	int start;
	int end;
	int count;
	int lseq;
}segmentInfo;

int countBamRec (bam1_t *b, bam_hdr_t *h, hts_itr_t *iter, htsFile *fp, int tid, int start, int end, int *l_qseq)
{
	bam1_core_t *c;
//	int id=omp_get_thread_num();
	int ret=0;
	int count=0;
	int lseq=0;
	int lastPos = 0;
	int firstPos = -1;
	while (ret >= 0)
	{
		c = &b->core;
		if (firstPos == -1)
			firstPos = c->pos;
		lastPos = c->pos;
//		if  (!((c->tid == tid) && (c->pos >= start) && (c->pos < end)))
//			printf ("c->tid %d c->pos %d (exp: tid %d start %d end %d)\n", c->tid, c->pos, tid, start, end);
		assert  (c->tid == tid);
		count ++;
		if (c->l_qseq > lseq) lseq = c->l_qseq;
		ret = sam_itr_next(fp, iter, b);
	}
//	printf ("(%d-%d): ", firstPos, lastPos);
	bam_hdr_destroy(h);
	bam_destroy1(b);
	*l_qseq = lseq;
	return count;
}
int main(int argc, char **argv)
{
	int i=0, j=0;
	if (argc != 4)
	{
		printf ("Usage: a.out <numThreads> <bamFile> <OutPrefix>\n");
		return 1;
	}
	int numThreads = atoi(argv[1]);
	char *bamFileName = argv[2];
	htsFile *fp1=hts_open(bamFileName, "r");
	bam_hdr_t *hdr1 = sam_hdr_read(fp1);
	
	int totalSegments=0;
  	int *numSegmentsPerContig = (int *)calloc (hdr1->n_targets, sizeof(int));
	assert (numSegmentsPerContig != NULL);
	for (i=0; i<hdr1->n_targets; i++)
	{
		numSegmentsPerContig[i] = hdr1->target_len[i]/SEGMENT_SIZE;
		if ((hdr1->target_len[i] % SEGMENT_SIZE) != 0) numSegmentsPerContig[i]++;
		totalSegments += numSegmentsPerContig[i];
	}
	printf ("totalSegments %d\n", totalSegments);
	int *contigWiseTotals = (int *)calloc (hdr1->n_targets, sizeof (int));
	assert (contigWiseTotals != NULL);

	segmentInfo *si = (segmentInfo *)malloc (totalSegments * sizeof (segmentInfo));
	int curSegNo = 0;
	for (i=0; i<hdr1->n_targets; i++)
	{
		for (j=0; j<numSegmentsPerContig[i]; j++)
		{
			si[curSegNo].contigNo = i;
			si[curSegNo].start = j*SEGMENT_SIZE;
			if (j==(numSegmentsPerContig[i]-1))
				si[curSegNo].end = hdr1->target_len[i];
			else
				si[curSegNo].end = (j+1)*SEGMENT_SIZE;
			si[curSegNo].count = 0;
			if (i==1)
				printf ("segNo %d start %d end %d\n", curSegNo, si[curSegNo].start, si[curSegNo].end);
			curSegNo++;
		}
	}
	#pragma omp parallel for private(i) num_threads(numThreads)
//	for (i=17; i<18/*hdr1->n_targets*/; i++)
	for (i=0; i<totalSegments; i++)
	{
		int tid = omp_get_thread_num();
		htsFile *fp=hts_open(bamFileName, "r");
		assert (fp!=NULL);
		hts_idx_t *idx = sam_index_load(fp, bamFileName);  // load the index
	        if (idx == NULL) {
	            printf("can't load index for \"%s\"", argv[1]);
	        }
		bam_hdr_t *hdr = sam_hdr_read(fp);
		bam1_t *b = (bam1_t *)calloc(1, sizeof(bam1_t));

		hts_itr_t *iter = sam_itr_queryi(idx, si[i].contigNo, si[i].start, si[i].end);
		if (iter == NULL) {
                	printf("can't parse region \"%s\"\n", hdr->target_name[i]);
            	}
		int ret = sam_itr_next(fp, iter, b);
		if (ret == -1) continue;
//		printf ("Tid %d Contig - %s \n",tid, hdr->target_name[i]);
		int count = countBamRec(b, hdr, iter, fp, si[i].contigNo, si[i].start, si[i].end, &(si[i].lseq));
		#pragma omp atomic
			contigWiseTotals[si[i].contigNo] += count;
		si[i].count = count;
		printf ("Tid %d Segment %d Contig %d (%s) NumReads %d, longest read %d\n", tid, i, si[i].contigNo, hdr1->target_name[si[i].contigNo], count, si[i].lseq);

		iter=NULL;
		hts_idx_destroy(idx);
		hts_close (fp);
	}
	for (i=0; i<hdr1->n_targets; i++)
	{
		printf ("%d %s %d\n", i, hdr1->target_name[i], contigWiseTotals[i]);
	}
	char fName[200];
	sprintf (fName, "%s.hist", bamFileName);
	FILE *bIdx = fopen (fName, "w");
	assert (bIdx != NULL);
	fwrite (si, sizeof(segmentInfo), totalSegments, bIdx);
	fclose (bIdx);

	hts_close (fp1);
	free (numSegmentsPerContig);
	free (si);
	free (contigWiseTotals);
	return 0;
}
