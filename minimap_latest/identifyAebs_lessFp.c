#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/stat.h>
#define NUMFP 1000

int numContigs = 0;
char **refFasta = NULL;
uint32_t BUFSIZE = 20000000u;
char *fastqBuffer  = NULL;
uint32_t buf_n=0, buf_m=0;

typedef struct
{
	int len;
	char name[100];
}contigInfo;

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

typedef struct 
{
        uint32_t pos;
        uint8_t seq[80];
        uint8_t quals[160];
        uint16_t flag;
        uint8_t qual;
        uint8_t  matchLen;
}fullRec;

typedef struct 
{
	int readId;
	char readName[100];
	int readLen;
	char seq[500];
	char qual[500];
	int qs;
	int qe;
	char strand;
	int contig;
	int contigLen;
	int rs;
	int re;

	int score;
	fullRec fR;
}pafRecord;

contigInfo *contigs=NULL;

void readFasta (char *fileName)
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
		printf ("%d - %s", i, contigs[i].name);
		while ((curInd < fSize) && (fastaContents[curInd]!='\n')) 
		{
//			printf ("%c", fastaContents[curInd]);
			curInd++;
		}
		printf ("\n");
		assert (curInd < fSize);
		curInd++;
		for (j=0; j<contigs[i].len; j++)
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
static int pafCompare (const void *p1, const void *q1)
{
	pafRecord * p = (pafRecord *)p1;
	pafRecord * q = (pafRecord *)q1;
	if ((p->score < q->score) || ((p->score == q->score) && (p->contig < q->contig)) || ((p->score == q->score) && (p->contig == q->contig) && (p->rs < q->rs)))
		return 1;
	else if ((p->score > q->score) || ((p->score == q->score) && (p->contig > q->contig)) || ((p->score == q->score) && (p->contig == q->contig) && (p->rs >  q->rs)))
		return -1;
	else
		return 0;
}
int calcMapQ (pafRecord *pr, int pr_n, int sub_n, int i)
{
	int mapQ=0;
	double identity = 1. - (double)(pr[0].readLen - pr[0].score) / 5 / pr[0].readLen;
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
void getNextRead (FILE *fpPaf, pafRecord *pr)
{
	if (fastqBuffer == NULL)	
	{
		fastqBuffer = (char *)malloc(BUFSIZE);
		assert (fastqBuffer != NULL);
		buf_m = fread (fastqBuffer, 1, BUFSIZE, fpPaf);
		printf ("After read1\n");
	}
	uint32_t i=0;
	for (i=buf_n; i<buf_m; i++)
	{
		if (fastqBuffer[i]=='\n') break;
	}
	if (i < buf_m)
	{
//		sscanf (&fastqBuffer[buf_n], "%d\t%s\t%d\t%s\t%s\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t", &pr->readId, pr->readName, &pr->readLen, pr->seq, pr->qual, &pr->qs, &pr->qe, &pr->strand, &pr->contig, &pr->contigLen, &pr->rs, &pr->re);
		pr->contig = 0;
		sscanf (&fastqBuffer[buf_n], "%d", &pr->readId);
		buf_n = i+1;
	}
	else if (!feof (fpPaf))
	{
		memmove (fastqBuffer, &fastqBuffer[buf_n], buf_m-buf_n);
		buf_m = buf_m-buf_n;
		buf_n = 0;
		size_t nread = fread (&fastqBuffer[buf_m], 1, BUFSIZE-buf_m, fpPaf);
		buf_m += nread;
		printf ("After read2 %lu\n", buf_m);
	}
	else
		pr->contig=-1;
	
}
void findAebRecords (char *pafFileName)
{
	pafRecord *pr1=NULL, *pr2=NULL;
	int pr1_m=0, pr2_m=0, pr1_n=0, pr2_n=0;

	int i=0, j=0;
	pafRecord pr;
	char line [1500];
	FILE *fp = fopen (pafFileName, "r");
	assert (fp != NULL);
	char fileName[200];
//	sprintf (fileName, "%s.aeb", pafFileName);
	FILE **fpAeb = (FILE **)calloc (NUMFP, sizeof (FILE *));
	assert (fpAeb != NULL);
	int *fpContigNo = (int *)malloc (NUMFP * sizeof(int));
	assert (fpContigNo != NULL);
	for (i=0; i<NUMFP; i++) fpContigNo[i]=-1;

	sprintf (fileName, "%s.rem.fq", pafFileName);
	FILE *fpOutFq = fopen (fileName, "w");
	assert (fpOutFq != NULL);

	fgets (line, 1500, fp);
//	getNextRead (fp, &pr);
	int minQual = 100, maxQual = 0;
	fullRec curFr, minFr;
	int curReadId, minReadId, curMatchCount, maxMatchCount, maxMatchCountOther, secondMaxMatchCount, sub_n;
	char minReadName[100], minReadSeq[200], minReadQual[200];
	char curReadName[100], curReadSeq[200], curReadQual[200];
	int firstRec=1;
#define DEBUG 0
	long minInsLen=100000, maxInsLen =0, numInsLen = 0, totInsLen=0, nextTarget=0;
	while (!feof (fp))
	{
		sscanf (line, "%d\t%s\t%d\t%s\t%s\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t", &pr.readId, pr.readName, &pr.readLen, pr.seq, pr.qual, &pr.qs, &pr.qe, &pr.strand, &pr.contig, &pr.contigLen, &pr.rs, &pr.re);
//		sscanf (line, "%d", &pr.readId );
		if ((pr.readId > nextTarget))
		{
			printf ("ReadId %d\n", pr.readId);
			nextTarget += 1000000;
		}
#if DEBUG
#else
		strcpy (curReadName, pr.readName);
		strcpy (curReadSeq, pr.seq);
		strcpy (curReadQual, pr.qual);
#endif
//		if (pr.readId>=80056) 
//			printf ("%d\n", pr.readId);
//			printf ("%d,%s,%d,%s,%s,%d,%d,%c,%d,%d,%d,%d \n ", pr.readId, pr.readName, pr.readLen, pr.seq, pr.qual, pr.qs, pr.qe, pr.strand, pr.contig, pr.contigLen, pr.rs, pr.re);
//		printf ("%s - ", line);
		int startPos=0;
		if (/*((pr.qe-pr.qs) == (pr.re-pr.rs)) &&*/ (pr.contigLen > (pr.rs-pr.qs+pr.readLen)))
		{
			float mapQ = 0;
			int numMatches=0;

			curReadId = pr.readId;
			curFr.matchLen = pr.readLen;
#if DEBUG
#else
			if (pr.strand == '+')
			{
				curFr.flag = 0;
				int rpos = pr.rs-pr.qs;
				startPos = rpos+1;
				curFr.pos = startPos;
				pr.rs=startPos;
//				if (pr.readId==80056)
//					printf ("REF -");
				for (i=0; i<pr.readLen; i++, rpos++)
				{
					if (pr.seq[i] > 96) pr.seq[i] -= 32;
					float qual = pr.qual[i] - 33;
					if ((pr.qual[i] - 33) > maxQual) maxQual = pr.qual[i] - 33;
					if ((pr.qual[i] - 33) < minQual) minQual = pr.qual[i] - 33;
//					qual = pow (10, -qual/10);
					if (pr.seq[i]==refFasta[pr.contig][rpos])
					{
						numMatches++;
					}
					else
					{
					}
					int seqInd = i/2;
					int ofs = i&1;
					if (ofs==0)
						curFr.seq[seqInd]=(uint8_t)seq_nt16_table[pr.seq[i]]<<4;
					else
						curFr.seq[seqInd] |= seq_nt16_table[pr.seq[i]];
					curFr.quals[i]=pr.qual[i] - 33;
//						curFr.seq[seqInd] |= (seq_nt16_table[pr.seq[i]]<<4);
				}
//				if (numMatches > 0)
//					printf ("%d,%s,%d,%s,%s,%d,%d,%c,%d,%d,%d,%d - NM:%d, MQ%.2f\n", pr.readId, pr.readName, pr.readLen, pr.seq, pr.qual, pr.qs, pr.qe, pr.strand, pr.contig, pr.contigLen, pr.rs, pr.re, numMatches);
			}
			else
			{
				curFr.flag = 16;
				int rpos = pr.rs - (pr.readLen - pr.qe);
				startPos = rpos+1;
				curFr.pos = startPos;
				pr.rs = startPos;
				for (i=(pr.readLen-1); i>=0; i--, rpos++)
				{
					if (pr.seq[i] > 96) pr.seq[i] -= 32;
					char rc = pr.seq[i];
					if (rc != 'N')
						rc = "TGAC"[(rc >> 1) & 3];
					float qual = pr.qual[i] - 33;
					if ((pr.qual[i] - 33) > maxQual) maxQual = pr.qual[i] - 33;
					if ((pr.qual[i] - 33) < minQual) minQual = pr.qual[i] - 33;
//					qual = pow (10, -qual/10);
					if (rc==refFasta[pr.contig][rpos])
					{
						numMatches++;
					}
					else
					{
					}
					int seqInd = (pr.readLen-1-i)/2;
					int ofs = (pr.readLen-1-i)&1;
					if (ofs==0)
						curFr.seq[seqInd]=(uint8_t)seq_nt16_table[rc]<<4;
					else
						curFr.seq[seqInd] |= seq_nt16_table[rc];
					curFr.quals[pr.readLen-1-i]=pr.qual[i] - 33;
				}
//					printf ("%s - NM:%d\n", line, numMatches);
			}
#endif
			pr.fR = curFr;
			pr.score = pr.readLen - 4*(pr.readLen-numMatches);
			if (pr.score < 0) pr.score = 0;
#if DEBUG
#else
			if ((curReadId != minReadId) || (firstRec==1))
			{
				if (curReadId % 2 == 0)
				{
					if (firstRec == 0)
					{//output aln for minReadId
//						printf ("pr1_n %d, pr2_n %d\n", pr1_n, pr2_n);
						if (pr1_n > 1)
							qsort (pr1, pr1_n, sizeof(pafRecord), pafCompare);
						if (pr2_n > 1)
							qsort (pr2, pr2_n, sizeof(pafRecord), pafCompare);

//						sub1_n=4; sub2_n=0;

						//Insert Length
						if ((pr1_n == 1) && (pr2_n == 1) && (pr1[0].contig == pr2[0].contig) && (pr1[0].score >= (pr1[0].readLen -20*4)) && (pr2[0].score >= (pr2[0].readLen -20*4)) && (abs (pr2[0].rs - pr1[0].rs) < 10000))
						{
							long insLen = (long)abs (pr2[0].rs - pr1[0].rs);
							if (insLen > maxInsLen) maxInsLen = insLen;
							if (insLen < minInsLen) minInsLen = insLen;
							totInsLen += insLen;
							numInsLen++;
						}
						if (pr1[0].readId % 1000000 == 0)
							printf ("ReadId %d Mean IS %.2lf, Min %ld, Max %ld, nObs %ld\n", pr1[0].readId, (double)totInsLen/numInsLen, minInsLen, maxInsLen, numInsLen);

						if (/*(pr1[0].contig == 10) && */(pr1[0].score >= (pr1[0].readLen -8*4)))
						{
							int posDiff=100000;
							i=1;
							while ((i<pr1_n) && (pr1[i].rs == pr1[0].rs)) i++;
							int sub1_n=0, sub2_n=0, j=0;
							//calc sub1_n
							for (j=0; j<pr1_n; j++)
							{
								if (j > i)
								{
									if ((pr1[j-1].rs != pr1[j].rs) && ((pr1[0].score - pr1[j].score) <= 7)) sub1_n++;
								}
//								if (pr1[].readId == 3102338)
//								printf ("%d: %d %c %d, %d, %d\n", j, pr1[j].readId, pr1[j].strand, pr1[j].contig, pr1[j].rs, pr1[j].score);
							}

							pr1[0].fR.qual = calcMapQ(pr1, pr1_n, sub1_n, i);
							pr2[0].fR.qual = 0;
							int i1=0;

							for (i=0; i<pr2_n; i++)
							{
//								pr1[0].fR.qual = mapQ_1;
								if ((pr2[i].contig == pr1[0].contig) && (abs (pr2[i].rs - pr1[0].rs) < posDiff) && (pr2[i].strand != pr1[0].strand) && (pr2[i].score >= (pr2[i].readLen -8*4)) && (posDiff > 300) && ((posDiff == 100000) || ((pr2[i1].score-pr2[i].score) <= 4*4)))
								{
//									if ((posDiff < 100000) && ((pr2[i1].score-pr2[i].score) > 3)) break;
/*									if (i==0)
										pr1[0].fR.qual = (mapQ_1 + mapQ_2)/2;
									else
										pr1[0].fR.qual = mapQ_1/2;
									pr2[i].fR.qual = pr1[0].fR.qual;*/
									i1=i;
									posDiff = abs (pr2[i].rs - pr1[0].rs);
//									if (posDiff < 300) break;
								}
							}
							if (posDiff < 100000)
							{
								if (i1==0)
								{
									int sub2Ind = 0, sub2PosDiff=100000;
									for (i=i1+1; i<pr2_n; i++)
									{
										if ((pr2[i-1].rs != pr2[i].rs) && (pr2[i].contig == pr1[0].contig) && (pr2[i].strand != pr1[0].strand) && (abs (abs (pr2[i].rs - pr1[0].rs) - abs(pr2[i1].rs - pr1[0].rs)) < sub2PosDiff) && ((sub2PosDiff == 100000) || (pr2[i].score >= (pr2[i].readLen -20*4))))
										{
											sub2Ind = i;
											sub2PosDiff = abs (abs (pr2[i].rs - pr1[0].rs) - abs(pr2[i1].rs - pr1[0].rs));
											if (pr2[i].score == pr2[i1].score) break;
										}
//										printf ("I %d Diff %d(%d-%d) sub2Ind %d sub2PosDiff %d\n", i, abs (pr2[i].rs - pr1[0].rs), pr2[i].rs, pr1[0].rs, sub2Ind, sub2PosDiff);
									}
									if (sub2Ind > 0) 
									{
//										printf ("sub2Ind %d sub2PosDiff %d\n", sub2Ind, sub2PosDiff);
									}
									else 
										sub2Ind = pr2_n;

									//calc sub2_n
									for (j=0; j<pr2_n; j++)
									{
										if (j >= sub2Ind)
										{
											if ((pr2[j-1].rs != pr2[j].rs) && ((pr2[0].score - pr2[j].score) <= 7)) sub2_n++;
										}	
//										if (pr1[].readId == 3102338)
//										printf ("%d: %d %c %d, %d, %d\n", j, pr2[j].readId, pr2[j].strand, pr2[j].contig, pr2[j].rs, pr2[j].score);
									}
									pr2[0].fR.qual = calcMapQ(pr2, pr2_n, sub2_n, sub2Ind);
								}
								i=i1;
//								printf ("pr1[0].contig %d\n", pr1[0].contig);
								if (fpContigNo[pr1[0].contig % NUMFP] != pr1[0].contig)
								{
									if (fpContigNo[pr1[0].contig % NUMFP] != -1)
									{
										fclose (fpAeb[pr1[0].contig % NUMFP]);
									}
									sprintf (fileName, "%s_%d.aeb", pafFileName, pr1[0].contig);
									fpAeb[pr1[0].contig % NUMFP] = fopen (fileName, "a");
									fpContigNo[pr1[0].contig % NUMFP] = pr1[0].contig;
								}
								assert (fpAeb[pr1[0].contig % NUMFP] != NULL);
								fwrite (&(pr1[0].fR), sizeof(fullRec), 1, fpAeb[pr1[0].contig % NUMFP]);
//								if (pr1[0].fR.qual > 0)
//								if (pr1[0].contig == 17)
								{
//									printf ("1-%d,%s/1,%d,%s,%s,%d,%d,%c,%d,%d,%d(%u),%d - NM:%d,MQ:%u\n", pr1[0].readId, pr1[0].readName, pr1[0].readLen, pr1[0].seq, pr1[0].qual, pr1[0].qs, pr1[0].qe, pr1[0].strand, pr1[0].contig, pr1[0].contigLen, pr1[0].rs, pr1[0].fR.pos, pr1[0].re, pr1[0].score, pr1[0].fR.qual);
								}
								if (fpContigNo[pr2[i].contig % NUMFP] != pr2[i].contig)
								{
									if (fpContigNo[pr2[i].contig % NUMFP] != -1)
									{
										fclose (fpAeb[pr2[i].contig % NUMFP]);
									}
									sprintf (fileName, "%s_%d.aeb", pafFileName, pr2[i].contig);
									fpAeb[pr2[i].contig % NUMFP] = fopen (fileName, "a");
									fpContigNo[pr2[i].contig % NUMFP] = pr2[i].contig;
								}
								assert (fpAeb[pr2[i].contig % NUMFP] != NULL);
								fwrite (&(pr2[i].fR), sizeof(fullRec), 1, fpAeb[pr2[i].contig % NUMFP]);
//								if (pr2[i].contig == 17)
//								if (pr2[i].fR.qual > 0)
								{
//									if (i!=0) printf ("*");
//									printf ("2-%d,%s/2,%d,%s,%s,%d,%d,%c,%d,%d,%d(%u),%d - NM:%d,MQ:%u\n", pr2[i].readId, pr2[i].readName, pr2[i].readLen, pr2[i].seq, pr2[i].qual, pr2[i].qs, pr2[i].qe, pr2[i].strand, pr2[i].contig, pr2[i].contigLen, pr2[i].rs, pr2[i].fR.pos, pr2[i].re, pr2[i].score, pr2[i].fR.qual);
								}
							}
							else
							{
								i=0;
//								pr2[i].fR.qual = mapQ_2;
								fprintf (fpOutFq, "@%s/1\n%s\n+\n%s\n@%s/2\n%s\n+\n%s\n", pr1[0].readName, pr1[0].seq, pr1[0].qual, pr2[0].readName, pr2[0].seq, pr2[0].qual);
//								fwrite (&(pr2[i].fR), sizeof(fullRec), 1, fpAeb);
//								if (mapQ_2 > 0)
//									printf ("3-%d,%s/2,%d,%s,%s,%d,%d,%c,%d,%d,%d,%d - NM:%d,MQ:%u,MQ1:%u,MQ2:%u\n", pr2[i].readId, pr2[i].readName, pr2[i].readLen, pr2[i].seq, pr2[i].qual, pr2[i].qs, pr2[i].qe, pr2[i].strand, pr2[i].contig, pr2[i].contigLen, pr2[i].rs, pr2[i].re, pr2[i].score, mapQ_2, mapQ_1, mapQ_2);
							}
						}
						else if (/*(pr2[0].contig == 10) &&*/ (pr2[0].score >= (pr2[0].readLen -8*4)))
						{
							int posDiff=100000;
							i=1;
							while ((i<pr2_n) && (pr2[i].rs == pr2[0].rs)) i++;
							int sub1_n=0, sub2_n=0, j=0;
							for (j=0; j<pr2_n; j++)
							{
								if (j > i)
								{
									if ((pr2[j-1].rs != pr2[j].rs) && ((pr2[0].score - pr2[j].score) <= 7)) sub2_n++;
								}
//								if (pr1[].readId == 3102338)
//								printf ("%d: %d %c %d, %d, %d\n", j, pr2[j].readId, pr2[j].strand, pr2[j].contig, pr2[j].rs, pr2[j].score);
							}
							pr2[0].fR.qual = calcMapQ(pr2, pr2_n, sub2_n, i);
							pr1[0].fR.qual = 0;
							int i1=0;
//							pr2[0].fR.qual = mapQ_2;
							for (i=0; i<pr1_n; i++)
							{
								if ((pr1[i].contig == pr2[0].contig) && (abs (pr1[i].rs - pr2[0].rs) < posDiff) && (pr1[i].strand != pr2[0].strand)  && (pr1[i].score >= (pr1[i].readLen -8*4)) && (posDiff > 300) && ((posDiff == 100000) || ((pr1[i1].score-pr1[i].score) <= 4*4)))
								{
//									if ((posDiff < 100000) && ((pr1[i1].score-pr2[i].score) > 3)) break;
/*									if (i==0)
										pr2[0].fR.qual = (mapQ_1 + mapQ_2)/2;
									else
										pr2[0].fR.qual = mapQ_2/2;
									pr1[i].fR.qual = pr2[0].fR.qual;*/
									i1=i;
									posDiff = abs (pr1[i].rs - pr2[0].rs);
//									if (posDiff < 300) break;
								}
							}
							if (posDiff < 100000)
							{
								if (i1==0)
								{
									int sub1Ind = 0, sub1PosDiff=100000;
									for (i=i1+1; i<pr1_n; i++)
									{
										if ((pr1[i-1].rs != pr1[i].rs) && (pr1[i].contig == pr2[0].contig) && (pr1[i].strand != pr2[0].strand) && (abs (abs (pr1[i].rs - pr2[0].rs) - abs(pr1[i1].rs - pr2[0].rs)) < sub1PosDiff) && ((sub1PosDiff == 100000) || (pr1[i].score >= (pr1[i].readLen -20*4))))
										{
											sub1Ind = i;
											sub1PosDiff = abs (abs (pr1[i].rs - pr2[0].rs) - abs(pr1[i1].rs - pr2[0].rs));
											if (pr1[i].score == pr1[i1].score) break;
										}
									}
									if (sub1Ind > 0) 
									{
//										printf ("sub1Ind %d sub1PosDiff %d\n", sub1Ind, sub1PosDiff);
									}
									else 
										sub1Ind = pr1_n;
									for (j=0; j<pr1_n; j++)
									{
										if (j >= sub1Ind)
										{
											if ((pr1[j-1].rs != pr1[j].rs) && ((pr1[0].score - pr1[j].score) <= 7)) sub1_n++;
										}
//										if (pr1[].readId == 3102338)
//										printf ("%d: %d %c %d, %d, %d\n", j, pr1[j].readId, pr1[j].strand, pr1[j].contig, pr1[j].rs, pr1[j].score);
									}
									pr1[0].fR.qual = calcMapQ(pr1, pr1_n, sub1_n, sub1Ind);
								}
								i=i1;
								if (fpContigNo[pr1[i].contig % NUMFP] != pr1[i].contig)
								{
									if (fpContigNo[pr1[i].contig % NUMFP] != -1)
									{
										fclose (fpAeb[pr1[i].contig % NUMFP]);
									}
									sprintf (fileName, "%s_%d.aeb", pafFileName, pr1[i].contig);
									fpAeb[pr1[i].contig % NUMFP] = fopen (fileName, "a");
									fpContigNo[pr1[i].contig % NUMFP] = pr1[i].contig;
								}
								assert (fpAeb[pr1[i].contig % NUMFP] != NULL);
								fwrite (&(pr1[i].fR), sizeof(fullRec), 1, fpAeb[pr1[i].contig % NUMFP]);
//								if (pr1[i].fR.qual > 0)
//								if (pr1[i].contig == 17)
								{
//									if (i!=0) printf ("*");
//									printf ("4-%d,%s/1,%d,%s,%s,%d,%d,%c,%d,%d,%d(%u),%d - NM:%d,MQ:%u\n", pr1[i].readId, pr1[i].readName, pr1[i].readLen, pr1[i].seq, pr1[i].qual, pr1[i].qs, pr1[i].qe, pr1[i].strand, pr1[i].contig, pr1[i].contigLen, pr1[i].rs, pr1[i].fR.pos, pr1[i].re, pr1[i].score, pr1[i].fR.qual);
								}
								if (fpContigNo[pr2[0].contig % NUMFP] != pr2[0].contig)
								{
									if (fpContigNo[pr2[0].contig % NUMFP] != -1)
										fclose (fpAeb[pr2[0].contig % NUMFP]);
									sprintf (fileName, "%s_%d.aeb", pafFileName, pr2[0].contig);
									fpAeb[pr2[0].contig % NUMFP] = fopen (fileName, "w");
									fpContigNo[pr2[0].contig % NUMFP] = pr2[0].contig;
								}
								assert (fpAeb[pr2[0].contig % NUMFP] != NULL);
								fwrite (&(pr2[0].fR), sizeof(fullRec), 1, fpAeb[pr2[0].contig % NUMFP]);
//								if (pr2[0].fR.qual > 0)
//								if (pr2[0].contig == 17)
								{
//									printf ("6-%d,%s/2,%d,%s,%s,%d,%d,%c,%d,%d,%d(%u),%d - NM:%d,MQ:%u\n", pr2[0].readId, pr2[0].readName, pr2[0].readLen, pr2[0].seq, pr2[0].qual, pr2[0].qs, pr2[0].qe, pr2[0].strand, pr2[0].contig, pr2[0].contigLen, pr2[0].rs, pr2[0].fR.pos, pr2[0].re, pr2[0].score, pr2[0].fR.qual);
								}
							}
							else
							{
								i=0;
//								pr1[i].fR.qual = mapQ_1;
								fprintf (fpOutFq, "@%s/1\n%s\n+\n%s\n@%s/2\n%s\n+\n%s\n", pr1[0].readName, pr1[0].seq, pr1[0].qual, pr2[0].readName, pr2[0].seq, pr2[0].qual);
//								fwrite (&(pr1[i].fR), sizeof(fullRec), 1, fpAeb);
/*								if (mapQ_1 > 0)
									printf ("5-%d,%s/1,%d,%s,%s,%d,%d,%c,%d,%d,%d,%d - NM:%d,MQ:%u,MQ1:%u,MQ2:%u\n", pr1[i].readId, pr1[i].readName, pr1[i].readLen, pr1[i].seq, pr1[i].qual, pr1[i].qs, pr1[i].qe, pr1[i].strand, pr1[i].contig, pr1[i].contigLen, pr1[i].rs, pr1[i].re, pr1[i].score, mapQ_1, mapQ_1, mapQ_2);*/
							}
						}
						else
						{
							fprintf (fpOutFq, "@%s/1\n%s\n+\n%s\n@%s/2\n%s\n+\n%s\n", pr1[0].readName, pr1[0].seq, pr1[0].qual, pr2[0].readName, pr2[0].seq, pr2[0].qual);
						}
					}
					firstRec=0;
					minReadId=curReadId;
					pr1_n=pr2_n=0;
				}
			}

			if (pr.readId % 2 == 0)
			{
				if (pr1_n == pr1_m)
				{
					pr1_m = (pr1_m == 0)? 16 : (pr1_m << 1);
//					printf ("pr1_m %d pr1_n %d pr1 %p\n", pr1_m, pr1_n, pr1);
					pr1 = (pafRecord *) realloc (pr1, pr1_m * sizeof (pafRecord));
					assert (pr1 != NULL);
				}
				memcpy (&pr1[pr1_n++], &pr, sizeof(pafRecord));
			}
			else
			{
				if (pr2_n == pr2_m)
				{
					pr2_m = (pr2_m == 0)? 16 : (pr2_m << 1);
					pr2 = (pafRecord *) realloc (pr2, pr2_m * sizeof (pafRecord));
					assert (pr2 != NULL);
				}
				memcpy (&pr2[pr2_n++], &pr, sizeof(pafRecord));
			}
#endif
			if (numMatches > 0)
			{
//				if ((pr.rs <= 185875) && (pr.re >= 185875) && (pr.contig==10)) 
//				if ((pr.readId == 429264069) || (pr.readId == 429264076) || (pr.readId == 559098375) || (pr.readId == 559098376) || (pr.readId == 559098441) || (pr.readId == 559098443) || (pr.readId == 646501105) || (pr.readId == 646501109) || (pr.readId == 646501110) || (pr.readId == 646501111) || (pr.readId == 646501182) || (pr.readId == 646501184) || (pr.readId == 646501193) || (pr.readId == 646501206) || (pr.readId == 646501216) || (pr.readId == 646501224) || (pr.readId == 646501225) || (pr.readId == 646501234) || (pr.readId == 706755002) || (pr.readId == 706755004))
//					printf ("%d,%s,%d,%s,%s,%d,%d,%c,%d,%d,%d,%d - NM:%d,MQ:%u\n", pr.readId, pr.readName, pr.readLen, pr.seq, pr.qual, pr.qs, pr.qe, pr.strand, pr.contig, pr.contigLen, pr.rs, pr.re, numMatches, curFr.qual);
			}
		}
		else
		{
//			fprintf (stderr/*fpOutFq*/, "@%s\n%s\n+\n%s\n", pr.readName, pr.seq, pr.qual);
		}
//		if (pr.readId==80056) 
//			printf ("\n");
		fgets (line, 1500, fp);
//		getNextRead (fp, &pr);
	}

	//output aln for minReadId
	
	fclose (fp);
	for (i=0; i<NUMFP;  i++)
	{
		if (fpAeb[i] != NULL) fclose (fpAeb[i]);
	}
	free (fpAeb);
	free (fpContigNo);
	fclose (fpOutFq);
	free (pr1);
	free (pr2);
}

int main(int argc, char **argv)
{
	int i=0;
	if (argc != 3)
	{
		printf ("Usage: a.out ref.fa aln.paf\n");
		return 1;
	}
	char *fastaFile = argv[1];
	char *pafFile = argv[2];

	char fileName[200];
	sprintf (fileName, "%s.fai", fastaFile);
	FILE *faiFp = fopen (fileName, "r");
	assert (faiFp != NULL);

	char line[200];
	fgets (line, 200, faiFp);
	while (!feof (faiFp))
	{
		numContigs ++;
		fgets (line, 200, faiFp);
	}
	printf ("numContigs %d\n", numContigs);
	fclose (faiFp);

	refFasta = (char **) malloc(numContigs * sizeof(char *));
	assert (refFasta != NULL);

	contigs = (contigInfo *)malloc (numContigs * sizeof (contigInfo));
	assert (contigs != NULL);
	faiFp = fopen (fileName, "r");
	assert (faiFp != NULL);
	for (i=0; i<numContigs; i++)
	{
		fgets (line, 200, faiFp);
		sscanf (line, "%s\t%d", contigs[i].name, &(contigs[i].len));
		refFasta[i] = (char *)malloc (contigs[i].len+1);
		assert (refFasta[i] != NULL);
	}

	readFasta (fastaFile);

	 int pid = getpid();
            char line1[2048];
        sprintf (line, "/proc/%d/status", pid);
            FILE *statFile = fopen(line, "r");
        assert (statFile != NULL);
            fgets (line, 2048, statFile);
        while (!feof (statFile))
            {
                    if (strstr(line,"VmPeak") || strstr(line,"VmHWM"))
                {
                        printf ("%s", line);
                    }
                fgets (line, 2048, statFile);
            }
            fclose (statFile);

	findAebRecords (pafFile);
	fclose (faiFp);
	free (contigs);
	for (i=0; i<numContigs; i++)
	{
		free (refFasta[i]);
	}
	free (refFasta);
	return 0;
}
