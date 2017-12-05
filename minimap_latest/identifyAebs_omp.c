#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/stat.h>
#include <omp.h>

int numContigs = 0;
char **refFasta = NULL;
int numThreads=1;
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
        uint8_t seq[60];
        uint8_t quals[120];
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
		printf ("%s", contigs[i].name);
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
	int p = ((pafRecord *)p1)->score;
	int q = ((pafRecord *)q1)->score;
	return (q-p);
}
void findAebRecords (char *pafFileName, int tid, size_t *offsets)
{
	pafRecord *pr1=NULL, *pr2=NULL;
	int pr1_m=0, pr2_m=0, pr1_n=0, pr2_n=0;
	int maxOutRecords = 200000;

	int i=0, j=0;
	pafRecord pr;
	char line [500];
	FILE *fp = fopen (pafFileName, "r");
	assert (fp != NULL);
	fseek (fp, offsets[tid], SEEK_SET);
	char fileName[200];
	sprintf (fileName, "%s_%d.aeb", pafFileName, tid);
	FILE *fpAeb = fopen (fileName, "w");
	assert (fpAeb != NULL);
	sprintf (fileName, "%s_%d.rem.fq", pafFileName, tid);
	FILE *fpOutFq = fopen (fileName, "w");
	assert (fpOutFq != NULL);

	fgets (line, 500, fp);
	int minQual = 100, maxQual = 0;
	fullRec curFr, minFr;
	fullRec *outRecs = (fullRec *)malloc(maxOutRecords * sizeof (fullRec));
	assert (outRecs != NULL);
	int curOutRec = 0;

	int curReadId, minReadId, curMatchCount, maxMatchCount, maxMatchCountOther, secondMaxMatchCount, sub_n;
	char minReadName[100], minReadSeq[200], minReadQual[200];
	char curReadName[100], curReadSeq[200], curReadQual[200];
	int firstRec=1;
	while (!feof (fp))
	{
		if (ftell(fp) >= offsets[tid+1]) break;
		sscanf (line, "%d\t%s\t%d\t%s\t%s\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t", &pr.readId, pr.readName, &pr.readLen, pr.seq, pr.qual, &pr.qs, &pr.qe, &pr.strand, &pr.contig, &pr.contigLen, &pr.rs, &pr.re);
		strcpy (curReadName, pr.readName);
		strcpy (curReadSeq, pr.seq);
		strcpy (curReadQual, pr.qual);
//		if (pr.readId>=80056) 
//			printf ("%d\n", pr.readId);
//			printf ("%d,%s,%d,%s,%s,%d,%d,%c,%d,%d,%d,%d - ", pr.readId, pr.readName, pr.readLen, pr.seq, pr.qual, pr.qs, pr.qe, pr.strand, pr.contig, pr.contigLen, pr.rs, pr.re);
//		printf ("%s - ", line);
		int startPos=0;
		if (/*((pr.qe-pr.qs) == (pr.re-pr.rs)) &&*/ (pr.contigLen > (pr.rs-pr.qs+pr.readLen)))
		{
			float mapQ = 0;
			int numMatches=0;

			curReadId = pr.readId;
			curFr.matchLen = pr.readLen;
//			memcpy (&(curFr.quals), &(pr.qual), pr.readLen);
			if (pr.strand == '+')
			{
				curFr.flag = 0;
				int rpos = pr.rs-pr.qs;
				startPos = rpos+1;
				curFr.pos = startPos;
				pr.rs=startPos;
/*				if (pr.readId>=80056)
					printf ("contig %d contigLen %d rpos %d-%d - ", pr.contig, pr.contigLen, rpos, rpos+pr.readLen);*/
//				if (pr.readId==80056)
//					printf ("REF -");
				for (i=0; i<pr.readLen; i++, rpos++)
				{
/*					if (pr.readId==80056)
					{
						printf ("%c", refFasta[pr.contig][rpos]);
					}*/
					if (pr.seq[i] > 96) pr.seq[i] -= 32;
					float qual = pr.qual[i] - 33;
					if ((pr.qual[i] - 33) > maxQual) maxQual = pr.qual[i] - 33;
					if ((pr.qual[i] - 33) < minQual) minQual = pr.qual[i] - 33;
//					qual = pow (10, -qual/10);
					if (pr.seq[i]==refFasta[pr.contig][rpos])
					{
/*						if (pr.readId == 7500758)
						{
							printf ("M:qual %c(%d) probability %f MapQ %f\n", pr.qual[i], pr.qual[i] - 33, qual, mapQ);
						}*/
						numMatches++;
					}
					else
					{
/*						if (pr.readId == 7500758)
						{
							printf ("N:qual %c(%d) probability %f MapQ %f\n", pr.qual[i], pr.qual[i] - 33, qual, mapQ);
						}*/
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
/*				if (pr.readId==80056) 
				{ 
					printf (" ");
				}*/
/*				if (numMatches > 0)
					printf ("%d,%s,%d,%s,%s,%d,%d,%c,%d,%d,%d,%d - NM:%d, MQ%.2f\n", pr.readId, pr.readName, pr.readLen, pr.seq, pr.qual, pr.qs, pr.qe, pr.strand, pr.contig, pr.contigLen, pr.rs, pr.re, numMatches);*/
			}
			else
			{
				curFr.flag = 16;
				int rpos = pr.rs - (pr.readLen - pr.qe);
				startPos = rpos+1;
				curFr.pos = startPos;
				pr.rs = startPos;
/*				if (pr.readId == 80067)
					printf ("rpos(%d)=%d(rs (%d) - [readLen(%d) - qe(%d)]\n", pr.readId, rpos, pr.rs, pr.readLen, pr.qe);*/
				for (i=(pr.readLen-1); i>=0; i--, rpos++)
				{
					if (pr.seq[i] > 96) pr.seq[i] -= 32;
					if (pr.seq[i] != 'N')
						pr.seq[i] = "TGAC"[(pr.seq[i] >> 1) & 3];
/*					if (pr.readId == 80067)
						printf ("r %c(%d) q %c(%d)\n", pr.seq[i], i, refFasta[pr.contig][rpos], rpos);*/
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
					int seqInd = (pr.readLen-1-i)/2;
					int ofs = (pr.readLen-1-i)&1;
					if (ofs==0)
						curFr.seq[seqInd]=(uint8_t)seq_nt16_table[pr.seq[i]]<<4;
					else
						curFr.seq[seqInd] |= seq_nt16_table[pr.seq[i]];
					curFr.quals[pr.readLen-1-i]=pr.qual[i] - 33;
				}
//					printf ("%s - NM:%d\n", line, numMatches);
			}
			pr.fR = curFr;
//			curFr.qual = (uint8_t)(-10*log10(mapQ));
			pr.score = numMatches;

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

						int sub1_n=0, sub2_n=0;
						for (i=0; i<pr1_n; i++)
						{
							if ((i>0) && ((pr1[0].score - pr1[i].score) > 5)) sub1_n++;
//							if (pr1[i].readId == 3102338)
								printf ("%d - %d, %d, %d\n", pr1[i].readId, pr1[i].contig, pr1[i].rs, pr1[i].score);
						}
						for (i=0; i<pr2_n; i++)
						{
							if ((i>0) && ((pr2[0].score - pr2[i].score) > 5)) sub2_n++;
//							if (pr1[i].readId == 3102338)
								printf ("%d - %d, %d, %d\n", pr2[i].readId, pr2[i].contig, pr2[i].rs, pr2[i].score);
						}
						int mapQ_1=0, mapQ_2=0;
						double identity = 1. - (double)(pr1[0].readLen - pr1[0].score) / 5 / pr1[0].readLen;
						double tmp = (double)3 / 4.61512;
						tmp *= identity * identity;
						if (pr1_n == 1)
							mapQ_1 = 60;
						else if (pr1[0].score == pr1[1].score)
							mapQ_1 = 0;
						else
						{
							int sub = (pr1[1].score < 19) ? 19: pr1[1].score;
							mapQ_1 = (int)(6.02 * (pr1[0].score - pr1[1].score) * tmp * tmp + .499);
							if (sub1_n > 0)
								mapQ_1 -= (int) (4.343 * log(sub1_n  + 1) + .499);
							if (mapQ_1 > 60) mapQ_1 = 60;
							if (mapQ_1 < 0) mapQ_1 =0;
						}
						identity = 1. - (double)(pr2[0].readLen - pr2[0].score) / 5 / pr2[0].readLen;
						tmp = (double)3 / 4.61512;
						tmp *= identity * identity;
						if (pr2_n == 1)
							mapQ_2 = 60;
						else if (pr2[0].score == pr2[1].score)
							mapQ_2 = 0;
						else
						{
							int sub = (pr2[1].score < 19) ? 19: pr2[1].score;
							mapQ_2 = (int)(6.02 * (pr2[0].score - pr2[1].score) * tmp * tmp + .499);
							if (sub2_n > 0)
								mapQ_2 -= (int) (4.343 * log(sub2_n  + 1) + .499);
							if (mapQ_2 > 60) mapQ_2 = 60;
							if (mapQ_2 < 0) mapQ_2 =0;
						}
//						printf ("mapQ_1 %d mapQ_2 %d sub1_n %d sub2_n %d\n", mapQ_1, mapQ_2, sub1_n, sub2_n);

						if ((pr1[0].contig == 10) && (pr1[0].score >= (pr1[0].readLen -4)))
						{
/*							for (i=0; i<pr2_n; i++)
							{
								pr1[0].fR.qual = mapQ_1;
								if ((pr2[i].contig == pr1[0].contig) && (abs (pr2[i].rs - pr1[0].rs) < 2000) && (pr2[i].strand != pr1[0].strand) && (pr2[i].score >= (pr2[i].readLen -10)))
								{
									if (i==0)
										pr1[0].fR.qual = (mapQ_1 + mapQ_2)/2;
									else
										pr1[0].fR.qual = mapQ_1/2;
									pr2[i].fR.qual = pr1[0].fR.qual;
									break;
								}
							}*/
							if (curOutRec == maxOutRecords)
							{
//								#pragma omp critical
								{
//									printf ("%d - curOutRec %d\n", tid, curOutRec);
									fwrite (outRecs, sizeof(fullRec), curOutRec, fpAeb);
								}
								curOutRec=0;
							}
							outRecs[curOutRec++] = pr1[0].fR;
//							if (pr1[0].fR.qual > 0)
							{
//								printf ("1-%d,%s/1,%d,%s,%s,%d,%d,%c,%d,%d,%d,%d - NM:%d,MQ:%u,MQ1:%u,MQ2:%u\n", pr1[0].readId, pr1[0].readName, pr1[0].readLen, pr1[0].seq, pr1[0].qual, pr1[0].qs, pr1[0].qe, pr1[0].strand, pr1[0].contig, pr1[0].contigLen, pr1[0].rs, pr1[0].re, pr1[0].score, pr1[0].fR.qual, mapQ_1, mapQ_2);
							}
/*							if (i < pr2_n)
							{
								if (curOutRec == maxOutRecords)
								{
//									#pragma omp critical
									{
//										printf ("%d - curOutRec %d\n", tid, curOutRec);
										fwrite (outRecs, sizeof(fullRec), curOutRec, fpAeb);
									}
									curOutRec=0;
								}
								outRecs[curOutRec++] = pr2[i].fR;
//								if (pr2[i].fR.qual > 0)
								{
//									if (i!=0) printf ("*");
//									printf ("2-%d,%s/2,%d,%s,%s,%d,%d,%c,%d,%d,%d,%d - NM:%d,MQ:%u,MQ1:%u,MQ2:%u\n", pr2[i].readId, pr2[i].readName, pr2[i].readLen, pr2[i].seq, pr2[i].qual, pr2[i].qs, pr2[i].qe, pr2[i].strand, pr2[i].contig, pr2[i].contigLen, pr2[i].rs, pr2[i].re, pr2[i].score, pr2[i].fR.qual, mapQ_1, mapQ_2);
								}
							}
							else
							{
								i=0;
								pr2[i].fR.qual = mapQ_2;
//								fwrite (&(pr2[i].fR), sizeof(fullRec), 1, fpAeb);
//								if (mapQ_2 > 0)
//									printf ("3-%d,%s/2,%d,%s,%s,%d,%d,%c,%d,%d,%d,%d - NM:%d,MQ:%u,MQ1:%u,MQ2:%u\n", pr2[i].readId, pr2[i].readName, pr2[i].readLen, pr2[i].seq, pr2[i].qual, pr2[i].qs, pr2[i].qe, pr2[i].strand, pr2[i].contig, pr2[i].contigLen, pr2[i].rs, pr2[i].re, pr2[i].score, mapQ_2, mapQ_1, mapQ_2);
							}*/
						}
						/*else */if ((pr2[0].contig == 10) && (pr2[0].score >= (pr2[0].readLen -4)))
						{
							pr2[0].fR.qual = mapQ_2;
/*							for (i=0; i<pr1_n; i++)
							{
								if ((pr1[i].contig == pr2[0].contig) && (abs (pr1[i].rs - pr2[0].rs) < 2000) && (pr1[i].strand != pr2[0].strand)  && (pr1[i].score >= (pr1[i].readLen -10)))
								{
									if (i==0)
										pr2[0].fR.qual = (mapQ_1 + mapQ_2)/2;
									else
										pr2[0].fR.qual = mapQ_2/2;
									pr1[i].fR.qual = pr2[0].fR.qual;
									break;
								}
							}
							if (i < pr1_n)
							{
								if (curOutRec == maxOutRecords)
								{
//									#pragma omp critical
									{
//										printf ("%d - curOutRec %d\n", tid, curOutRec);
										fwrite (outRecs, sizeof(fullRec), curOutRec, fpAeb);
									}
									curOutRec=0;
								}
								outRecs[curOutRec++] = pr1[i].fR;
//								if (pr1[i].fR.qual > 0)
								{
//									if (i!=0) printf ("*");
//									printf ("4-%d,%s/1,%d,%s,%s,%d,%d,%c,%d,%d,%d,%d - NM:%d,MQ:%u,MQ1:%u,MQ2:%u\n", pr1[i].readId, pr1[i].readName, pr1[i].readLen, pr1[i].seq, pr1[i].qual, pr1[i].qs, pr1[i].qe, pr1[i].strand, pr1[i].contig, pr1[i].contigLen, pr1[i].rs, pr1[i].re, pr1[i].score, pr1[i].fR.qual, mapQ_1, mapQ_2);
								}
							}
							else
							{
								i=0;
								pr1[i].fR.qual = mapQ_1;
//								fwrite (&(pr1[i].fR), sizeof(fullRec), 1, fpAeb);
//								if (mapQ_1 > 0)
//									printf ("5-%d,%s/1,%d,%s,%s,%d,%d,%c,%d,%d,%d,%d - NM:%d,MQ:%u,MQ1:%u,MQ2:%u\n", pr1[i].readId, pr1[i].readName, pr1[i].readLen, pr1[i].seq, pr1[i].qual, pr1[i].qs, pr1[i].qe, pr1[i].strand, pr1[i].contig, pr1[i].contigLen, pr1[i].rs, pr1[i].re, pr1[i].score, mapQ_1, mapQ_1, mapQ_2);
							}*/
							if (curOutRec == maxOutRecords)
							{
//								#pragma omp critical
								{
//									printf ("%d - curOutRec %d\n", tid, curOutRec);
									fwrite (outRecs, sizeof(fullRec), curOutRec, fpAeb);
								}
								curOutRec=0;
							}
							outRecs[curOutRec++] = pr2[0].fR;
//							if (pr2[0].fR.qual > 0)
							{
//								printf ("6-%d,%s/2,%d,%s,%s,%d,%d,%c,%d,%d,%d,%d - NM:%d,MQ:%u,MQ1:%u,MQ2:%u\n", pr2[0].readId, pr2[0].readName, pr2[0].readLen, pr2[0].seq, pr2[0].qual, pr2[0].qs, pr2[0].qe, pr2[0].strand, pr2[0].contig, pr2[0].contigLen, pr2[0].rs, pr2[0].re, pr2[0].score, pr2[0].fR.qual, mapQ_1, mapQ_2);
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


			if (numMatches > 0)
			{
//				if ((pr.rs <= 185875) && (pr.re >= 185875) && (pr.contig==10)) 
//				if ((pr.readId == 429264069) || (pr.readId == 429264076) || (pr.readId == 559098375) || (pr.readId == 559098376) || (pr.readId == 559098441) || (pr.readId == 559098443) || (pr.readId == 646501105) || (pr.readId == 646501109) || (pr.readId == 646501110) || (pr.readId == 646501111) || (pr.readId == 646501182) || (pr.readId == 646501184) || (pr.readId == 646501193) || (pr.readId == 646501206) || (pr.readId == 646501216) || (pr.readId == 646501224) || (pr.readId == 646501225) || (pr.readId == 646501234) || (pr.readId == 706755002) || (pr.readId == 706755004))
//					printf ("%d,%s,%d,%s,%s,%d,%d,%c,%d,%d,%d,%d - NM:%d,MQ:%u\n", pr.readId, pr.readName, pr.readLen, pr.seq, pr.qual, pr.qs, pr.qe, pr.strand, pr.contig, pr.contigLen, pr.rs, pr.re, numMatches, curFr.qual);
			}
		}
		else
		{
//			fprintf (fpOutFq, "@%s\n%s\n+\n%s\n", pr.readName, pr.seq, pr.qual);
		}
//		if (pr.readId==80056) 
//			printf ("\n");
		fgets (line, 500, fp);
	}
	if (curOutRec > 0)
	{
		#pragma omp critical
		{ 
//			printf ("%d - curOutRec %d\n", tid, curOutRec);
			fwrite (outRecs, sizeof(fullRec), curOutRec, fpAeb);
		}
	}
	free (outRecs);
	//output aln for minReadId
//	printf ("pr1_n %d, pr2_n %d\n", pr1_n, pr2_n);
	int sub1_n=0, sub2_n=0;
	if (pr1_n > 1)
		qsort (pr1, pr1_n, sizeof(pafRecord), pafCompare);
	if (pr2_n > 1)
		qsort (pr2, pr2_n, sizeof(pafRecord), pafCompare);
	for (i=0; i<pr1_n; i++)
	{
		if ((i>0) && ((pr1[0].score - pr1[i].score) > 5)) sub1_n++;
//		printf ("%d - %d, %d, %d\n", pr1[i].readId, pr1[i].contig, pr1[i].rs, pr1[i].score);
	}
	for (i=0; i<pr2_n; i++)
	{
		if ((i>0) && ((pr2[0].score - pr2[i].score) > 5)) sub2_n++;
//		printf ("%d - %d, %d, %d\n", pr2[i].readId, pr2[i].contig, pr2[i].rs, pr2[i].score);
	}
	int mapQ_1=0, mapQ_2=0;
	double identity = 1. - (double)(pr1[0].readLen - pr1[0].score) / 5 / pr1[0].readLen;
	double tmp = (double)3 / 4.61512;
	tmp *= identity * identity;
	if (pr1_n == 1)
		mapQ_1 = 60;
	else
	{
		int sub = (pr1[1].score < 19) ? 19: pr1[1].score;
		mapQ_1 = (int)(6.02 * (pr1[0].score - pr1[1].score) * tmp * tmp + .499);
		if (sub1_n > 0)
			mapQ_1 -= (int) (4.343 * log(sub1_n  + 1) + .499);
		if (mapQ_1 > 60) mapQ_1 = 60;
		if (mapQ_1 < 0) mapQ_1 =0;
	}
	identity = 1. - (double)(pr2[0].readLen - pr2[0].score) / 5 / pr2[0].readLen;
	tmp = (double)3 / 4.61512;
	tmp *= identity * identity;
	if (pr2_n == 1)
		mapQ_2 = 60;
	else
	{
		int sub = (pr2[1].score < 19) ? 19: pr2[1].score;
		mapQ_2 = (int)(6.02 * (pr2[0].score - pr2[1].score) * tmp * tmp + .499);
		if (sub2_n > 0)
			mapQ_2 -= (int) (4.343 * log(sub2_n  + 1) + .499);
		if (mapQ_2 > 60) mapQ_2 = 60;
		if (mapQ_2 < 0) mapQ_2 =0;
	}
//	printf ("mapQ_1 %d mapQ_2 %d sub1_n %d, sub2_n %d\n", mapQ_1, mapQ_2, sub1_n, sub2_n);

	
	fclose (fp);
	fclose (fpAeb);
	fclose (fpOutFq);
	free (pr1);
	free (pr2);
}

int main(int argc, char **argv)
{
	int i=0;
	if (argc != 4)
	{
		printf ("Usage: a.out <numThreads> ref.fa aln.paf\n");
		return 1;
	}
	numThreads = atoi(argv[1]);
	char *fastaFile = argv[2];
	char *pafFile = argv[3];

	char fileName[200];
	sprintf (fileName, "%s.fai", fastaFile);
	FILE *faiFp = fopen (fileName, "r");
	assert (faiFp != NULL);

	char line[500];
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

	size_t *offsets = (size_t *)calloc (numThreads+1, sizeof(size_t));
	struct stat buf;
	int res = stat (pafFile, &buf);
	assert (res==0);
	size_t fSize = buf.st_size;
	offsets[numThreads] = fSize;
	FILE *fp = fopen (pafFile, "r");
	pafRecord pr1, pr2;
	for (i=1; i<numThreads; i++)
	{
		size_t ofs = fSize*i/numThreads;
		fseek (fp, ofs, SEEK_SET);
		fgets (line, 500, fp);
		fgets (line, 500, fp);
		sscanf (line, "%d\t%s\t%d\t%s\t%s\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t", &pr1.readId, pr1.readName, &pr1.readLen, pr1.seq, pr1.qual, &pr1.qs, &pr1.qe, &pr1.strand, &pr1.contig, &pr1.contigLen, &pr1.rs, &pr1.re);
		
		ofs = ftell (fp);
		fgets (line, 500, fp);
		sscanf (line, "%d\t%s\t%d\t%s\t%s\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t", &pr2.readId, pr2.readName, &pr2.readLen, pr2.seq, pr2.qual, &pr2.qs, &pr2.qe, &pr2.strand, &pr2.contig, &pr2.contigLen, &pr2.rs, &pr2.re);
		while (!feof(fp) && ((pr2.readId == pr1.readId) || (pr2.readId & 1 == 1)))
		{
//			printf ("while %d\n", pr2.readId);
			pr1 = pr2;
			ofs = ftell (fp);
			fgets (line, 500, fp);
			if (!feof(fp))
				sscanf (line, "%d\t%s\t%d\t%s\t%s\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t", &pr2.readId, pr2.readName, &pr2.readLen, pr2.seq, pr2.qual, &pr2.qs, &pr2.qe, &pr2.strand, &pr2.contig, &pr2.contigLen, &pr2.rs, &pr2.re);
		}
		assert (!feof(fp));
		offsets[i] = ofs;
		printf ("offsets[%d]=%lu\n", i, ofs);
	}
	fclose (fp);
	readFasta (fastaFile);
/*	sprintf (fileName, "%s.aeb", pafFile);
	FILE *fpAeb = fopen (fileName, "w");
	assert (fpAeb != NULL);
	fclose (fpAeb);*/

	#pragma omp parallel num_threads (numThreads)
	{
		findAebRecords (pafFile, omp_get_thread_num(), offsets);
	}
	fclose (faiFp);
	free (contigs);
	for (i=0; i<numContigs; i++)
	{
		free (refFasta[i]);
	}
	free (refFasta);
	free (offsets);
	return 0;
}
