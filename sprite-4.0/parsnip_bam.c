#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <errno.h>
#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef USE_OMP
#include <omp.h>
#endif
#include <math.h>
#include <getopt.h>
#include <stdint.h>
#include <sys/stat.h>
#include <assert.h>
#include "sam.h"
#include "hts.h"
#include "bgzf.h"

char sam_prefix[500];
int NUM_CONTIGS=93;
char fai_file[200];
char vcf_file[200];
char fileName[200];
char fasta_file[200];

int nextPos=-1;
long posList[333216];
int chunkSize = 262144*16;
int numThreads = 1;
FILE *fp_fasta;
int totalAlnFiles = 0;
int numtasks=1, rank=0;
int WINDOW_SIZE=50000;
#define SEGMENT_SIZE 1000000
int maxBufferRecords=1000000;
typedef struct {
	int contigNo;
	int start;
	int end;
	int count;
}segmentInfo;
segmentInfo *si=NULL;
//FILTERS
/*
 * MQ - (avg_mapq >= 20)
 * DP - (pos_table[ind].count[4] > 1)
 * SB - (r1>.1)
 * AAF - (maxPer > 20)
 * AAC - (max_per > 1)
 */ 
int MQ=20, DP=1, AAF=20, AAC=1;
float SB=.1;

char **ref_chars=NULL;
typedef struct ref_pos
{
   uint32_t count[5];   
   uint32_t readPos[4];
   uint16_t fwdCount;//#alt alleles in fwd reads
   uint16_t fwdTotal;//#total fwd reads
   uint16_t revCount;
   uint16_t revTotal;
   uint32_t mapQ;
}ref_pos;
typedef struct 
{
	int contigNo;
	char *fileName;
	long numRecords;
	long seqLength;
	long numAeb;
	long numAib;
        long REF_START;
}fileInfo;
fileInfo *fileList=NULL;

typedef struct
{
        uint32_t pos;
        uint8_t *seq;
	uint8_t *quals;
	uint16_t flag;
        uint8_t qual;
        uint8_t  matchLen;
}fullRec;

typedef struct
{
        uint32_t pos;
	uint16_t flag;
        uint8_t *seq;
	uint8_t *quals;
        uint8_t qual;
        uint8_t n_cigar;
        int *cigar;
}otherRec;

typedef struct
{
 	fullRec *fRecs;
	otherRec *oRecs;
	int numFullRecs;
	int numOtherRecs;
}window;

long maxSeqLen=0;
int startContig=-1, endContig=-1;
int *procAssign=NULL;
int *startPos=NULL;
int *contigPosCount=NULL;

static int sizecmpfn(const void *a, const void *b) {
    fileInfo *av = (fileInfo *) a;
    fileInfo *bv = (fileInfo *) b;
    return av->numRecords-bv->numRecords;
}
void readFasta (int contigNo, FILE *fp_fasta, int len)
{
	long ind=0;
	fseek(fp_fasta, fileList[contigNo].REF_START, SEEK_SET);
	int lenToRead = len + (len/50)+100;
/*	if (contigStr==NULL)
	{
		printf ("Alloc fastaIp Str %ld\n", maxSeqLen + (maxSeqLen/50)+110);
	        contigStr = (char *)malloc(maxSeqLen + (maxSeqLen/50)+110);
		assert (contigStr != NULL);
	}*/
	char *contigStr = (char *)malloc(lenToRead+10);
        assert (contigStr != NULL);
        fread (contigStr, 1, lenToRead, fp_fasta);
	fclose (fp_fasta);
	int curIndex=0;
        for (ind=0; ind<len; ind++)
        {
	    char ch=contigStr[curIndex++];
            if (ch=='\n')
            {
                ind--;
                continue;
            }
            if (islower(ch))
                ch=toupper(ch);
            ref_chars[contigNo][ind]=ch;
        }
	free (contigStr);
}
inline int processOther (otherRec *rec, int winStart, int contigNo, ref_pos *posTable, int *indelPosition)
{	
	int ret=0;
	int i=0,j=0;
	int pos=0;
	int startPos = rec->pos-1;
	startPos = startPos - winStart;
	int ind=0;
	char Seq[110];
//	if (rec->n_cigar > 10) rec->n_cigar=10;
	for (i=0; i<rec->n_cigar; i++)
	{
		int cigarType=rec->cigar[i] & 0xF;
		int len = rec->cigar[i]>>4;
		switch(cigarType)
		{
			case 0://Match
				for(j=0;j<len;j++)
				{
					if (rec->qual >0)
						posTable[startPos].count[4]++;
					int curInd=pos/2;
				        int subInd= pos&1;
					ind=(rec->seq[curInd]>>((1-subInd)*4))&15;
//                               		ind=(rec->seq[curInd]>>((3-subInd)*2))&3;
					Seq[j]="=ACMGRSVTWYHKDBN"[ind];
					ind = (Seq[j]>>1) & 3;
					if ((Seq[j] == 'A') || (Seq[j] == 'C') || (Seq[j] == 'G') || (Seq[j] == 'T'))
					{	
//						if ((startPos+winStart) == 1200014) printf ("%d-%c(SAMAlnPos %d)\n", j, Seq[j], rec->pos);
						if (rec->qual >0)
							posTable[startPos].count[ind]++;
						if (rec->flag & 0x10)
							posTable[startPos].revTotal++;
						else
							posTable[startPos].fwdTotal++;
						posTable[startPos].readPos[ind]+=pos;
						if (ref_chars[contigNo][startPos+winStart] != Seq[j])
						{
							posTable[startPos].mapQ+=rec->qual;
							if (rec->flag & 0x10)
								posTable[startPos].revCount++;
							else
								posTable[startPos].fwdCount++;
						}
					}
					pos++;
					startPos++;
				}
				Seq[j]=0;
				break;
			case 1://Insert
				ret=1;
/*				if (*numIndels==*indelArraySize)
				{
					*indelArraySize=(*indelArraySize==0)?1000:*indelArraySize<<1;
					*indelPositions = (int *)realloc (*indelPositions, *indelArraySize*sizeof(int));
					assert (*indelPositions != NULL);
				}
				(*indelPositions)[*numIndels]=startPos+winStart;
				*numIndels = *numIndels + 1;*/
				*indelPosition=startPos+winStart;
				pos+=len;
				break;
			case 2://Delete
				ret=1;
/*				if (*numIndels==*indelArraySize)
				{
					*indelArraySize=(*indelArraySize==0)?1000:*indelArraySize<<1;
					*indelPositions = (int *)realloc (*indelPositions, *indelArraySize*sizeof(int));
					assert (*indelPositions != NULL);
				}
				(*indelPositions)[*numIndels]=startPos+winStart;
				*numIndels = *numIndels + 1;*/
				*indelPosition=startPos+winStart;
				for(j=0;j<len;j++)
				{
					if (rec->qual >0)
						posTable[startPos].count[4]++;
					startPos++;
				}
				break;
			case 4://Soft
				pos+=len;
				break;
			case 5://H
//				pos+=len;
				break;
		}
	}
	return ret;
}
void outputVcfRecs(int startPos, FILE *fp_vcf, ref_pos *pos_table, char *cName, int seqLen, FILE *fn, int contigNo){
	long ind=0;
//	fprintf (stderr, "Start Window: %d-%d\n", startPos, startPos + seqLen-1);
	for (ind=0; ind<seqLen; ind++)
        {
        	char ch=ref_chars[contigNo][ind+startPos];
                uint16_t max_per=0, ref_per=0;
		int pref,palt;
                char max_char=' ';
		switch(ch)
                {
                case 'A':
                    ref_per=pos_table[ind].count[0];
                    if (pos_table[ind].count[1] > pos_table[ind].count[3])
                    {
                       max_char='C';
                       max_per=pos_table[ind].count[1];
		       palt=(pos_table[ind].count[1]==0)?0:pos_table[ind].readPos[1]/pos_table[ind].count[1];
                    }
                    else
                    {
                       max_char='G';
                       max_per=pos_table[ind].count[3];
		       palt=(pos_table[ind].count[3]==0)?0:pos_table[ind].readPos[3]/pos_table[ind].count[3];
                    }
                    if (pos_table[ind].count[2] > max_per)
                    {
                       max_char='T';
                       max_per=pos_table[ind].count[2];
		       palt=(pos_table[ind].count[2]==0)?0:pos_table[ind].readPos[2]/pos_table[ind].count[2];
                    }
		    pref=(pos_table[ind].count[0]==0)?0:pos_table[ind].readPos[0]/pos_table[ind].count[0];
                    break;
		case 'C':
                    ref_per=pos_table[ind].count[1];
                    if (pos_table[ind].count[0] > pos_table[ind].count[3])
                    {
                       max_char='A';
                       max_per=pos_table[ind].count[0];
		       palt=(pos_table[ind].count[0]==0)?0:pos_table[ind].readPos[0]/pos_table[ind].count[0];
                    }
                    else
                    {
                       max_char='G';
                       max_per=pos_table[ind].count[3];
		       palt=(pos_table[ind].count[3]==0)?0:pos_table[ind].readPos[3]/pos_table[ind].count[3];
                    }
                    if (pos_table[ind].count[2] > max_per)
                    {
                       max_char='T';
                       max_per=pos_table[ind].count[2];
		       palt=(pos_table[ind].count[2]==0)?0:pos_table[ind].readPos[2]/pos_table[ind].count[2];
                    }
                    break;
		    pref=(pos_table[ind].count[1]==0)?0:pos_table[ind].readPos[1]/pos_table[ind].count[1];
                 case 'T':
                    ref_per=pos_table[ind].count[2];
                    if (pos_table[ind].count[0] > pos_table[ind].count[1])
                    {
                       max_char='A';
                       max_per=pos_table[ind].count[0];
		       palt=(pos_table[ind].count[0]==0)?0:pos_table[ind].readPos[0]/pos_table[ind].count[0];
                    }
                    else
                    {
                       max_char='C';
                       max_per=pos_table[ind].count[1];
		       palt=(pos_table[ind].count[1]==0)?0:pos_table[ind].readPos[1]/pos_table[ind].count[1];
                    }
                    if (pos_table[ind].count[3] > max_per)
                    {
                       max_char='G';
                       max_per=pos_table[ind].count[3];
		       palt=(pos_table[ind].count[3]==0)?0:pos_table[ind].readPos[3]/pos_table[ind].count[3];
                    }
		    pref=(pos_table[ind].count[2]==0)?0:pos_table[ind].readPos[2]/pos_table[ind].count[2];
                    break;
		case 'G':
                    ref_per=pos_table[ind].count[3];
                    if (pos_table[ind].count[0] > pos_table[ind].count[1])
                    {
                       max_char='A';
                       max_per=pos_table[ind].count[0];
		       palt=(pos_table[ind].count[0]==0)?0:pos_table[ind].readPos[0]/pos_table[ind].count[0];
                    }
                    else
                    {
                       max_char='C';
                       max_per=pos_table[ind].count[1];
		       palt=(pos_table[ind].count[1]==0)?0:pos_table[ind].readPos[1]/pos_table[ind].count[1];
                    }
                    if (pos_table[ind].count[2] > max_per)
                    {
                       max_char='T';
                       max_per=pos_table[ind].count[2];
		       palt=(pos_table[ind].count[2]==0)?0:pos_table[ind].readPos[2]/pos_table[ind].count[2];
                    }
		    pref=(pos_table[ind].count[3]==0)?0:pos_table[ind].readPos[3]/pos_table[ind].count[3];
                    break;
		}
		float maxPer = (float)max_per*100/pos_table[ind].count[4];
//		double maxPer1 = 1-((double)max_per/pos_table[ind].count[4]);
//		double qualPhred = -10 * (log(maxPer1) / log(10));
		int gt1=0;
//                if (((startPos+ind+1) == 1241903) || ((startPos+ind+1) == 4132573) || ((startPos+ind+1) == 4158057) || ((startPos+ind+1) == 6607766))
//                {
//			int avg_mapq = 0;
//			if (max_per > ref_per)
//				avg_mapq = (pos_table[ind].count[4]==0)?0:(pos_table[ind].mapQ/(pos_table[ind].count[4]-ref_per));
//	                printf("%s\t%ld\t.\t%c\t%c\t%d\t.", cName, startPos+ind+1, ch, max_char, avg_mapq);
//        	        printf("\tDP=%d;ACTG=%d,%d,%d,%d;FC=%u;RC=%u;FT=%u;RT=%u;AAF=%.2f\n", pos_table[ind].count[4], pos_table[ind].count[0], pos_table[ind].count[1], pos_table[ind].count[2], pos_table[ind].count[3], pos_table[ind].fwdCount, pos_table[ind].revCount, pos_table[ind].fwdTotal, pos_table[ind].revTotal, maxPer);
//                }
                if ((maxPer > AAF) && (max_char != ch) && ((max_per > AAC)) && (pos_table[ind].count[4] > DP))
                {
		   if (maxPer > 80) gt1=1;
//				   float fwdRatio = (pos_table[ind].fwdTotal>0)?((float)(pos_table[ind].fwdCount)/(pos_table[ind].fwdTotal)):0;
//				   float revRatio = (pos_table[ind].revTotal>0)?((float)(pos_table[ind].revCount)/(pos_table[ind].revTotal)):0;
				   int isZero = 0;
				   if ((pos_table[ind].fwdTotal ==0) || (pos_table[ind].revTotal ==0)) isZero=1;
//				   float ratioOfRatio = (fwdRatio<revRatio)?(fwdRatio/revRatio):(revRatio/fwdRatio);
				   float r1 = (pos_table[ind].fwdCount < pos_table[ind].revCount)?((float)pos_table[ind].fwdCount/pos_table[ind].revCount):((float)pos_table[ind].revCount/pos_table[ind].fwdCount);
				   float r2 = (pos_table[ind].fwdCount < pos_table[ind].revCount)?((float)pos_table[ind].fwdTotal/pos_table[ind].revTotal):((float)pos_table[ind].revTotal/pos_table[ind].fwdTotal);
				int sbLowStrandCount = (pos_table[ind].fwdCount < pos_table[ind].revCount)?pos_table[ind].fwdTotal:pos_table[ind].revTotal;
				int avg_mapq = (pos_table[ind].count[4]==0)?0:(pos_table[ind].mapQ/(pos_table[ind].count[4]-ref_per));
				if ( ((r1>SB) || ((r1==0) && (isZero==1))) && (avg_mapq >= MQ))
				{	
		                   fprintf(fp_vcf, "%s\t%ld\t.\t%c\t%c\t%d\t.", cName, startPos+ind+1, ch, max_char, avg_mapq);
		                   fprintf(fp_vcf, "\tDP=%d;ACTG=%d,%d,%d,%d;FC=%u;RC=%u;FT=%u;RT=%u;SB=%.2f;SBT=%.3f;SBLow=%d;AAF=%.2f", pos_table[ind].count[4], pos_table[ind].count[0], pos_table[ind].count[1], pos_table[ind].count[2], pos_table[ind].count[3], pos_table[ind].fwdCount, pos_table[ind].revCount, pos_table[ind].fwdTotal, pos_table[ind].revTotal, r1, r2, sbLowStrandCount, maxPer);
/*				   if (maxPer < 25) fprintf(fp_vcf, "20");
				   else if (maxPer <30) fprintf(fp_vcf, "25");
				   else if (maxPer <35) fprintf(fp_vcf, "30");
				   else fprintf(fp_vcf, "35");*/
		                   fprintf(fp_vcf, "\tGT\t%d/1\n", gt1);
				}
/*				
			}*/
                }
        }
}
void parsnip_init(int seqLen, char *fasta_file, int contig)
{
/*        posTable = (ref_pos *)malloc(seqLen * sizeof(ref_pos));
        assert(posTable != NULL);
        bzero(posTable, seqLen * sizeof(ref_pos));*/
        sprintf (fileName, "%s_c%d", vcf_file, fileList[contig].contigNo);
        fp_fasta = fopen(fasta_file, "r");
	assert (fp_fasta != NULL);
        readFasta (contig, fp_fasta, seqLen);
	assert ((WINDOW_SIZE%10)==0);
}
void parsnip_aeb(fullRec *Fr, int winStart, int contigNo, ref_pos *posTable)
{
        int l;
        char Seq[150];
//	Fr->qual *= ;
	float q = (float)Fr->qual;
//	q = q*60/41;
//	Fr->qual = (uint8_t) q;
//	if (Fr->qual > 0) Fr->qual = 50;
        totalAlnFiles=1;
                long pos=Fr->pos-1;
		pos = pos - winStart;
                for (l=0; l<Fr->matchLen; l++)
                {
                        int curInd=l/2;
                        int subInd= l&1;
                        uint8_t ind=(Fr->seq[curInd]>>((1-subInd)*4))&15;
			Seq[l]="=ACMGRSVTWYHKDBN"[ind];
			ind = (Seq[l]>>1) & 3;
			if ((Seq[l] == 'A') || (Seq[l] == 'C') || (Seq[l] == 'G') || (Seq[l] == 'T'))
			{
//			Seq[l]="ACTG"[ind];
				if (Fr->qual > 0)
				{
		                        posTable[pos].count[ind]++;
//				if ((winStart+pos) == 1200014) printf ("F%d-%c SAMPos %d WinStart %d pos %d\n", l, Seq[l], Fr->pos, winStart, pos);
        		                posTable[pos].count[4]++;
				}
				if (Fr->flag & 0x10)
        	                        posTable[pos].revTotal++;
                	        else
                        	        posTable[pos].fwdTotal++;
				if (Seq[l] != ref_chars[contigNo][pos+winStart])
        	                {
                	                posTable[pos].mapQ+=Fr->qual;
					if (pos+winStart == 4158056)
						printf ("Pos 4158056 MQ %d\n", Fr->qual);
                        	        if (Fr->flag & 0x10)
                                	     	posTable[pos].revCount++;
	                                else
        	                                posTable[pos].fwdCount++;
                	        }
			}
			pos++;
		}
		Seq[l]=0;
}			
int parsnip_aib(otherRec *Or, int winStart, int contigNo, ref_pos *posTable, int *indelPosition)
{
        totalAlnFiles=1;
	return processOther(Or, winStart, contigNo, posTable, indelPosition);
}
void parsnip_write(char *cName, int seqLen, int startPos, FILE *fp_vcf, int contigNo, ref_pos *posTable)
{
        if (totalAlnFiles > 0)
        {
                assert(fp_vcf!=NULL);
                outputVcfRecs(startPos, fp_vcf, posTable, cName, seqLen, NULL, contigNo);
        }
}

int getSeqLen (otherRec *oR)
{
	int i=0, seqLen=0;
	for (i=0; i<oR->n_cigar; i++)
	{
		uint32_t len = oR->cigar[i]>>4;
		seqLen += len;
	}
	return seqLen;
}

void processSegment (int segNo, int rank, int tid, FILE *fp_vcf, char *bamFile)
{
	ref_pos *posTable=NULL;
//	printf ("[%d]process segment %d Contig %d:%d-%d Count %d\n", rank, segNo, si[segNo].contigNo, si[segNo].start, si[segNo].end, si[segNo].count);
	posTable = (ref_pos *)malloc((SEGMENT_SIZE+SEGMENT_SIZE/10) * sizeof(ref_pos));
//	printf ("Alloc postable for contig %d - %ld entries\n", si[segNo].contigNo, WINDOW_SIZE+WINDOW_SIZE/10/*fileList[i].seqLength*/);
	assert(posTable!=NULL);
	bzero(posTable, (SEGMENT_SIZE+SEGMENT_SIZE/10)/*fileList[i].seqLength*/ * sizeof(ref_pos));

	int indelPosition=0, curIndelPos=0;

	bam1_core_t *c;
	fullRec FR;
	otherRec OR;
	char *data;

	size_t max_k = 0, k = 0;
	bam1_t **buf=NULL;

	htsFile *fp=hts_open(sam_prefix, "r");
	assert (fp!=NULL);
	hts_idx_t *idx = sam_index_load(fp, sam_prefix);
	if (idx == NULL) {
		printf("can't load index for \"%s\"", sam_prefix);
	}
	bam_hdr_t *hdr = sam_hdr_read(fp);
	bam1_t *b = (bam1_t *)calloc(1, sizeof(bam1_t));
//	printf ("get iter: segNo %d contigNo %d, start %d, end %d\n", segNo, si[segNo].contigNo, si[segNo].start, si[segNo].end);
	hts_itr_t *iter = sam_itr_queryi(idx, si[segNo].contigNo, si[segNo].start, si[segNo].end);
	if (iter == NULL) {
		printf("can't parse region \"%s\"\n", hdr->target_name[si[segNo].contigNo]);
	}
	int ret = sam_itr_next(fp, iter, b);
	int count=0;
	int winStart=-1, winEnd=-1;
	while (ret >= 0)
	{
		c = &b->core;
		data = (char *)b->data;
		int *cigar=(int *)(data+c->l_qname);
		uint8_t *seq = bam_get_seq(b);
		uint8_t *qualStr = bam_get_qual(b);
		if (winStart==-1) winStart = c->pos;
		winEnd = c->pos;
		if ((c->n_cigar==1) && (cigar[0]%16==0))
		{
			//AEB
			FR.pos = c->pos+1;
			FR.qual = c->qual;
			FR.flag = c->flag;
			FR.matchLen = cigar[0]/16;
			FR.quals = qualStr;
			FR.seq = seq;
			assert (c->l_qseq == FR.matchLen);
			parsnip_aeb (&FR, winStart, si[segNo].contigNo, posTable);
		}
		else
		{
			OR.pos = c->pos+1;
			OR.flag = c->flag;
			OR.qual = c->qual;
			OR.n_cigar = c->n_cigar;
			OR.quals = qualStr;
			OR.seq = seq;
			OR.cigar = cigar;
			int ret = parsnip_aib(&OR, winStart, si[segNo].contigNo, posTable, &indelPosition);
			if (ret)
				curIndelPos = indelPosition;
		}
		if ((c->pos <= curIndelPos) && ((c->pos + c->l_qseq + 100) > curIndelPos) && (c->pos > si[segNo].start))
		{
			if (k == max_k) {
				size_t old_max = max_k;
				max_k = max_k? max_k<<1 : 0x10000;
				buf = (bam1_t**)realloc(buf, max_k * sizeof(bam1_t*));
				memset(buf + old_max, 0, sizeof(bam1_t*) * (max_k - old_max));
			}
			buf[k] = (bam1_t*)calloc(1, sizeof(bam1_t));
			bam1_t *b1 = buf[k];
			memcpy (b1, b, sizeof(bam1_t));
			b1->data = (uint8_t*)malloc (b->m_data);
			assert (b1->data != NULL);
			memcpy (b1->data, b->data, b->l_data);	
			k++;
		}
		count++;
		ret = sam_itr_next(fp, iter, b);
	}
	if (count > 0)
	{
		printf ("[%d, %d]segNo %d Contig %d Count %d ExpCount %d winStart %d winEnd %d(%d) NumIndels %d\n", rank, tid, segNo, si[segNo].contigNo, count, si[segNo].count, winStart, winEnd, si[segNo].end, k);
		if (winStart < si[segNo].start)
		{
			bzero (posTable, (si[segNo].start-winStart)*sizeof(ref_pos));
		}
		parsnip_write(hdr->target_name[si[segNo].contigNo], si[segNo].end-winStart, winStart, fp_vcf, si[segNo].contigNo, posTable);
		int i;
		samFile* fpOut;
		fpOut = sam_open(bamFile, "wb1");
		if (fpOut == NULL) return;
		sam_hdr_write(fpOut, hdr);
		hts_set_threads(fpOut, 1);
		for (i = 0; i < k; ++i)
			sam_write1(fpOut, hdr, buf[i]);
		sam_close(fpOut);
	}
	else if (segNo == rank)
	{
/*		samFile* fpOut;
		fpOut = sam_open(bamFile, "wb1");
		if (fpOut == NULL) return;
		printf ("segNo %d write header\n", segNo);
		sam_hdr_write(fpOut, hdr);
		sam_close(fpOut);*/
	}
	for (k = 0; k < max_k; ++k) {
		if (!buf[k]) continue;
		free(buf[k]->data);
		free(buf[k]);
	}
	free(buf);
	hts_idx_destroy (idx);
	hts_close (fp);
	hts_itr_destroy (iter);
	bam_hdr_destroy(hdr);
	bam_destroy1(b);
	
	free (posTable);
}

int main(int argc, char **argv)
{
        struct stat buf;
        int i=0,j=0;
	static struct option long_options[] =
	{
		{"MQ", required_argument,       0, 'm'},
		{"DP", required_argument,       0, 'd'},
		{"SB", required_argument,       0, 's'},
		{"AAF", required_argument,       0, 'f'},
		{"AAC", required_argument,       0, 'c'},
		{"nthreads", required_argument,       0, 't'},
		{0, 0, 0, 0}
	};
	int c=0;
	while (1)
        {
		int option_index = 0;
		c = getopt_long (argc, argv, "m:s:d:f:c:t:", long_options, &option_index);
		if (c==-1) break;
		switch (c)
		{
			case 'm':
				MQ=atoi(optarg);
				break;
			case 'd':
				DP=atoi(optarg);
				break;
			case 's':
				SB=(float)atof(optarg);
				break;
			case 'f':
				AAF=(float)atoi(optarg);
				break;
			case 'c':
				AAC=(float)atoi(optarg);
				break;
			case 't':
				numThreads=(int)atoi(optarg);
				break;
		}
	}
        if (optind + 1 >= argc || optind + 3 != argc) 
        {
		fprintf (stderr, "\n");
                fprintf (stderr, "Usage: parsnip [options] <path to fasta file> <sorted BAM file> <output vcf file prefix>\n");
                fprintf (stderr, "Options:\n\n");
		fprintf (stderr, "	-t INT        number of threads [1]\n");
		fprintf (stderr, "	--MQ=INT      minimum average SNP mapping quality filter [20]\n");
		fprintf (stderr, "	--DP=INT      minimum read depth filter [2]\n");
		fprintf (stderr, "	--SB=REAL     minimum strand bias filter [.1]\n");
		fprintf (stderr, "	--AAF=INT     minimum alternate allele frequency filter [20]\n");
		fprintf (stderr, "	--AAC=INT     minimum alternate allele count filter [2]\n");
                return 0;
        }
        strcpy(fasta_file, argv[optind]);
        sprintf (fai_file, "%s.fai", fasta_file);
	FILE *fp_fai = fopen(fai_file, "r");
        strcpy (sam_prefix, argv[optind+1]);
        strcpy (vcf_file, argv[optind+2]);

#ifdef USE_MPI
        MPI_Init(&argc,&argv);
        MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
	printf ("rank %d\n",rank);
        char line[200];
        long length=0,start=0;

        NUM_CONTIGS=0;
        fgets(line, 200, fp_fai);
        while (!feof(fp_fai))
        {
                NUM_CONTIGS++;
                fgets(line, 200, fp_fai);
        }
        fclose(fp_fai);
        fileList = (fileInfo *)malloc(NUM_CONTIGS * sizeof(fileInfo));
        fp_fai = fopen(fai_file, "r");
        fgets(line, 200, fp_fai);
        int contigNo=0;
        char contigName[100];
        char fileName[500];
        long totalRecords=0;

//        sprintf(fileName, "%s.hist", sam_prefix);
//	int status = stat(fileName, &buf);
//	assert (status == 0);
//	int numSegments = buf.st_size/sizeof(segmentInfo);
//	si = (segmentInfo *)malloc (numSegments * sizeof(segmentInfo));
//	assert (si != NULL);
//	FILE *fpSeg = fopen (fileName, "r");
//	fread (si, sizeof(segmentInfo), numSegments, fpSeg);
//	fclose (fpSeg);
	
        while (!feof(fp_fai))
        {
                sscanf(line, "%s\t%ld\t%ld", contigName, &length, &start);
		fileList[contigNo].contigNo=contigNo;
                int fLen = strlen(contigName);
                fileList[contigNo].fileName = (char *)malloc(fLen+1);
                strcpy(fileList[contigNo].fileName, contigName);
                fileList[contigNo].numRecords = 0;
                fileList[contigNo].seqLength = length;
                fileList[contigNo].REF_START = start;
                if (maxSeqLen < length)
                        maxSeqLen = length;
                contigNo++;
                fgets(line, 200, fp_fai);
        }
        fclose(fp_fai);
	ref_chars = (char **)calloc (NUM_CONTIGS, sizeof(char *));
	printf ("[%d]reading FASTA.. \n", rank);
	
        int totalSegments=0;
        int *numSegmentsPerContig = (int *)calloc (NUM_CONTIGS, sizeof(int));
        assert (numSegmentsPerContig != NULL);

	for (i=0; i<NUM_CONTIGS; i++)
	{
		ref_chars[i] = (char *) malloc (fileList[i].seqLength * sizeof (char));
		assert (ref_chars[i] != NULL);
		parsnip_init(fileList[i].seqLength,fasta_file, i);
                numSegmentsPerContig[i] = fileList[i].seqLength/SEGMENT_SIZE;
                if ((fileList[i].seqLength % SEGMENT_SIZE) != 0) numSegmentsPerContig[i]++;
                totalSegments += numSegmentsPerContig[i];
	}
	printf ("[%d]Done \n totalSegments %d\n", rank, totalSegments);
	si = (segmentInfo *)malloc (totalSegments * sizeof (segmentInfo));
        int curSegNo = 0;
        for (i=0; i<NUM_CONTIGS; i++)
        {
                for (j=0; j<numSegmentsPerContig[i]; j++)
                {
                        si[curSegNo].contigNo = i;
                        si[curSegNo].start = j*SEGMENT_SIZE;
                        if (j==(numSegmentsPerContig[i]-1))
                                si[curSegNo].end = fileList[i].seqLength;
                        else
                                si[curSegNo].end = (j+1)*SEGMENT_SIZE;
                        si[curSegNo].count = 0;
                        if (i==1)
                                printf ("segNo %d start %d end %d\n", curSegNo, si[curSegNo].start, si[curSegNo].end);
                        curSegNo++;
                }
        }
	free (numSegmentsPerContig);

//        qsort(fileList, NUM_CONTIGS, sizeof(fileInfo), sizecmpfn);
        long recordsPerProc = totalRecords / numtasks;
//        procAssign = (int *)malloc(sizeof(int) * (numtasks+1));
//        bzero(procAssign, sizeof(int) * (numtasks+1));
//        procAssign[numtasks]=numSegments;
        long prefixSum=si[0].count;
        int segInd=0;
        long index=0;

	int totalSize=0;
//	for (i=0; i<numSegments; i++)
//		totalSize += si[i].count;
	recordsPerProc = totalSegments/numtasks;
	
	int tid=0;
      	struct timeval st;
        gettimeofday(&st, NULL);
#ifdef USE_OMP
	#pragma omp parallel for num_threads (numThreads) schedule(dynamic)
#endif
		for (i=rank; i<totalSegments; i+=numtasks)
		{
#ifdef USE_OMP
			tid = omp_get_thread_num();
#endif
//			if (si[i].contigNo<=1)
			{
			      	char suffix[500];
        			sprintf (suffix, "%s_C%d_%d.vcf", vcf_file, si[i].contigNo, i);
				FILE *fp_vcf=fopen (suffix, "w");
				assert (fp_vcf != NULL);
        			sprintf (suffix, "%s_C%d_%d_indelReads.bam", vcf_file, si[i].contigNo, i);
				processSegment (i, rank, tid, fp_vcf, suffix);
				fclose (fp_vcf);
			}
		}
	char suffix[500];
	char *vcfBuf=NULL;
	size_t vcfSize=0;
#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
       	struct timeval et;
        gettimeofday(&et, NULL);
	if (rank==0)
	{
		long elapsed = (et.tv_sec * 1000000 + et.tv_usec) - (st.tv_sec * 1000000 + st.tv_usec);
		printf ("VarCall completed in %.2lf sec\nCombining VCF files...\n", (double)elapsed/1000000);
		sprintf (suffix, "%s.vcf", vcf_file);
		FILE *fp_vcf = fopen(suffix, "w");
		assert (fp_vcf != NULL);
                fprintf (fp_vcf, "##fileformat=VCFv4.1\n##source=SPRITE3.0\n##reference=%s\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth at the locus\">\n##INFO=<ID=ACTG,Number=4,Type=Integer,Description=\"Count for each allele type\">\n##INFO=<ID=FC,Number=1,Type=Integer,Description=\"Alternate allele count in forward strands\">\n##INFO=<ID=RC,Number=1,Type=Integer,Description=\"Alternate allele count in reverse strands\">\n##INFO=<ID=FT,Number=1,Type=Integer,Description=\"Total forward strands\">\n##INFO=<ID=RT,Number=1,Type=Integer,Description=\"Total reverse strands\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tvcf\n",fasta_file);
		int nread=0;
		for (i=0; i<totalSegments; i++)
		{
			sprintf (suffix, "%s_C%d_%d.vcf", vcf_file, si[i].contigNo, i);
			int ret = stat(suffix, &buf);
			if (ret == 0)
			{
				if (buf.st_size > vcfSize)
				{
					vcfSize = buf.st_size;
					vcfBuf = (char *)realloc (vcfBuf, vcfSize);
					assert (vcfBuf != NULL);
				}
				FILE *tmpVcfFp = fopen (suffix, "r");
				nread = fread(vcfBuf, 1, buf.st_size, tmpVcfFp);
				assert (nread == buf.st_size);
				fwrite (vcfBuf, 1, nread, fp_vcf);
				fclose (tmpVcfFp);
				int stat=unlink(suffix);
				if (stat==-1)
					printf ("Cannot remove temporary VCF file %s. Error code %d\n", suffix, errno);
			}
		}
		fclose (fp_vcf);
		printf ("Done.\n");
		printf ("Generating BAM file with INDEL reads...\n");
		sprintf (suffix, "%s_indelReads.bam", vcf_file);
		samFile* fpOut;
		fpOut = sam_open(suffix, "wb1");
		if (fpOut == NULL) return 1;
		hts_set_threads(fpOut, 1);
		int first=1;
		for (i=0; i<totalSegments; i++)
		{
			sprintf (suffix, "%s_C%d_%d_indelReads.bam", vcf_file, si[i].contigNo, i);
			int ret = stat(suffix, &buf);
			if (ret == 0)
			{
				BGZF *fp=bgzf_open(suffix, "r");
				assert (fp!=NULL);
				bam_hdr_t *hdr = bam_hdr_read(fp);
				bam1_t *b = (bam1_t *)calloc(1, sizeof(bam1_t));
				printf ("Writing %s\n", suffix);
				if (first)
				{
					first=0;
					sam_hdr_write(fpOut, hdr);
				}
/*				if (buf.st_size > vcfSize)
				{
					vcfSize = buf.st_size;
					vcfBuf = (char *)realloc (vcfBuf, vcfSize);
					assert (vcfBuf != NULL);
				}*/
//				FILE *tmpVcfFp = fopen (suffix, "r");
//				nread = fread(vcfBuf, 1, buf.st_size, tmpVcfFp);
//				assert (nread == buf.st_size);
//				fwrite (vcfBuf, 1, nread, fp_vcf);
//				fclose (tmpVcfFp);
				while (bam_read1(fp, b) >= 0) 
				{
					sam_write1(fpOut, hdr, b);
				}
				bgzf_close (fp);
				bam_hdr_destroy(hdr);
				bam_destroy1(b);
//				int stat=unlink(suffix);
//				if (stat==-1)
//					printf ("Cannot remove temporary VCF file %s. Error code %d\n", suffix, errno);
			}
		}
//		fclose (fp_vcf);
		sam_close(fpOut);
		printf ("Done.\n");
		if (vcfBuf != NULL)
			free (vcfBuf);
	}
/*	for (i=1; i<1; i++)
	{
		printf ("Contig %d\n", i);
        	struct timeval tm;
	        gettimeofday(&tm, NULL);


		int initDone=0;
		totalAlnFiles=0;
		sprintf(fileName, "%s_%s_sorted.aeb", sam_prefix, fileList[i].fileName);
		FILE *fp_aeb = fopen(fileName, "r");
		sprintf(fileName, "%s_%s_sorted.aib", sam_prefix, fileList[i].fileName);
		FILE *fp_aib = fopen(fileName, "r");

		int maxAebRecs = 0, curAebRec = 0, maxAibRecs = 0, curAibRec = 0, eofAeb=0, eofAib=0;
		fullRec *fullBuffer = NULL;
		otherRec *otherBuffer = NULL;

		

		fullRec *Fr = NULL;//(fullRec *)malloc(sizeof(fullRec));
//        	assert(Fr != NULL);
       		otherRec *Or = NULL;//(otherRec *)malloc(sizeof(otherRec));
//        	assert(Or != NULL);
		int curWindow=0, winStart=0, winEnd=WINDOW_SIZE, indelWindow=0, indelWindowCount=0, maxWindowSize=0;
		if (fp_aeb != NULL)
		{
			assert (fileList[i].numAeb > 0);
			if (fileList[i].numAeb < maxBufferRecords)
				maxAebRecs = fileList[i].numAeb;
			else
				maxAebRecs = maxBufferRecords;

			fullBuffer = (fullRec *)malloc (maxAebRecs * sizeof (fullRec));
			assert (fullBuffer != NULL);
			
			assert (fread (fullBuffer, sizeof(fullRec), maxAebRecs, fp_aeb) == maxAebRecs);
			if ((maxAebRecs < maxBufferRecords) ||  (feof(fp_aeb)))
				eofAeb=1;
//			winStart = Fr->pos-1;
//			winEnd = winStart + Fr->matchLen -2;
			initDone = 1;
		}
		else
		{
			printf ("File %s not found\n", fileName);
		}
		if (fp_aib != NULL)
                {
			assert (fileList[i].numAib > 0);
			if (fileList[i].numAib  < maxBufferRecords)
				maxAibRecs = fileList[i].numAib;
			else
				maxAibRecs = maxBufferRecords;
			otherBuffer = (otherRec *) malloc (maxAibRecs * sizeof (otherRec));
			assert (otherBuffer != NULL);

			assert (fread (otherBuffer, sizeof(otherRec), maxAibRecs, fp_aib) == maxAibRecs);
			if ((maxAibRecs < maxBufferRecords) || (feof(fp_aib)))
				eofAib=1;
//			int seqLen = getSeqLen(Or);
			if (initDone==0)
			{
				parsnip_init(fileList[i].seqLength,fasta_file, i);
				initDone=1;
			}
		}
		else
		{
			printf ("File %s not found\n", fileName);
		}
		printf ("WinStart %d WinEnd %d\n", winStart, winEnd);
 		while ((fp_aeb && (curAebRec < maxAebRecs)) || (fp_aib && (curAibRec < maxAibRecs)))
		{
			if (fp_aeb && (curAebRec < maxAebRecs))
			{
				Fr = &fullBuffer[curAebRec];
				if ((fp_aib && (curAibRec < maxAibRecs) && (Fr->pos <= otherBuffer[curAibRec].pos)) || (!fp_aib) || (fp_aib && (curAibRec >= maxAibRecs)))
				{
					while (Fr->pos > winEnd)
					{
//						if (!indelWindow)
						{//writeWindow
			                        	parsnip_write(fileList[i].fileName, winEnd-winStart, winStart, fp_vcf, i);
						}
//TODO::bzero
						bzero (posTable, WINDOW_SIZE*sizeof(ref_pos));
						memmove (posTable, &posTable[WINDOW_SIZE], (WINDOW_SIZE/10)*sizeof(ref_pos));
						bzero (&posTable[WINDOW_SIZE], (WINDOW_SIZE/10)*sizeof(ref_pos));
						winStart = winEnd;
						winEnd = winStart + WINDOW_SIZE;
						curWindow++;
						if (indelWindow)
							indelWindowCount++;
						indelWindow=0;
					}
					if (Fr->pos > winEnd)
					{
						printf ("Fr->pos %d winEnd %d\n", Fr->pos, winEnd);
						assert (0);
					}
//					if (Fr->qual > 0)
					{
						parsnip_aeb (Fr, winStart, i);
					}
//					if ((Fr->pos <= 42354966) && ((Fr->pos + Fr->matchLen-1) >= 42354966))
//						printf ("F: %d-%d, Win(%d-%d)\n", Fr->pos, Fr->pos + Fr->matchLen-1, winStart, winEnd);
//					else
//						return 1;
					curAebRec++;
					if (curAebRec == maxAebRecs)
					{
						if (!eofAeb)
						{
							maxAebRecs=fread (fullBuffer, sizeof(fullRec), maxAebRecs, fp_aeb);
							if (feof(fp_aeb))
								eofAeb=1;
							curAebRec=0;
						}
					}
				}
				else
				{
					Or = &otherBuffer[curAibRec];
					int seqLen = getSeqLen(Or);
					while (Or->pos > winEnd)
					{
//						if (!indelWindow)
						{//writeWindow
                        				parsnip_write(fileList[i].fileName, winEnd-winStart, winStart, fp_vcf, i);
						}
						bzero (posTable, WINDOW_SIZE*sizeof(ref_pos));
						memmove (posTable, &posTable[WINDOW_SIZE], (WINDOW_SIZE/10)*sizeof(ref_pos));
						bzero (&posTable[WINDOW_SIZE], (WINDOW_SIZE/10)*sizeof(ref_pos));
						winStart = winEnd;
						winEnd = winStart + WINDOW_SIZE;
						if (indelWindow)
							indelWindowCount++;
						indelWindow=0;
						curWindow++;
					}
					if (Or->pos > winEnd)
					{
						printf ("Or->pos %d winEnd %d\n", Or->pos, winEnd);
						assert (0);
					}
//					if (Or->qual > 0)
						indelWindow |= parsnip_aib (Or, winStart, i);
//					if ((Or->pos <= 42354966) && ((Or->pos + getSeqLen(Or) -1) >= 42354966))
//						printf ("V: %d-%d, Win(%d-%d)\n", Or->pos, Or->pos + getSeqLen(Or) -1, winStart, winEnd);
//					else
//						return 1;
					curAibRec++;
					if (curAibRec == maxAibRecs)
					{
						if (!eofAib)
						{
							maxAibRecs=fread (otherBuffer, sizeof(otherRec), maxAibRecs, fp_aib);
							if (feof(fp_aib))
								eofAib=1;
							curAibRec=0;
						}
					}
				}
			}
			else if (fp_aib && (curAibRec < maxAibRecs))
			{

				Or = &otherBuffer[curAibRec];
				int seqLen = getSeqLen(Or);
				while (Or->pos > winEnd)
				{
//					if (!indelWindow)
					{//writeWindow
						parsnip_write(fileList[i].fileName, winEnd-winStart, winStart, fp_vcf, i);
					}
					bzero (posTable, WINDOW_SIZE*sizeof(ref_pos));
					memmove (posTable, &posTable[WINDOW_SIZE], (WINDOW_SIZE/10)*sizeof(ref_pos));
					bzero (&posTable[WINDOW_SIZE], (WINDOW_SIZE/10)*sizeof(ref_pos));
					winStart = winEnd;
					winEnd = winStart + WINDOW_SIZE;
					if (indelWindow)
						indelWindowCount++;
					indelWindow=0;
					curWindow++;
				}
				if (Or->pos > winEnd)
				{
					printf ("Or->pos %d winEnd %d\n", Or->pos, winEnd);
					assert (0);
				}
//				if (Or->qual > 0)
					indelWindow |= parsnip_aib (Or, winStart, i);
//				if ((Or->pos <= 42354966) && ((Or->pos + getSeqLen(Or) -1) >= 42354966))
//					printf ("V: %d-%d, Win(%d-%d)\n", Or->pos, Or->pos + getSeqLen(Or) -1, winStart, winEnd);
//				else
//					return 1;
				curAibRec++;
				if (curAibRec == maxAibRecs)
				{
					if (!eofAib)
					{
						maxAibRecs=fread (otherBuffer, sizeof(otherRec), maxAibRecs, fp_aib);
						if (feof(fp_aib))
							eofAib=1;
						curAibRec=0;
					}
				}
			}
		}
		printf ("outOfWhile aeb %d/%d eof %d, aib %d/%d eof %d\n", curAebRec, maxAebRecs, eofAeb, curAibRec, maxAibRecs, eofAib);
//                parsnip_write(fileList[i].fileName, fileList[i].seqLength, 0, fp_vcf);
                parsnip_write(fileList[i].fileName, winEnd-winStart, winStart, fp_vcf, i);
		bzero (posTable, WINDOW_SIZE*sizeof(ref_pos));
		memmove (posTable, &posTable[WINDOW_SIZE], (WINDOW_SIZE/10)*sizeof(ref_pos));
		bzero (&posTable[WINDOW_SIZE], (WINDOW_SIZE/10)*sizeof(ref_pos));
		winStart = winEnd;
		winEnd = winStart + WINDOW_SIZE;
                parsnip_write(fileList[i].fileName, winEnd-winStart, winStart, fp_vcf, i);
		if (fp_aeb) fclose (fp_aeb);
		if (fp_aib) fclose (fp_aib);
		if (fullBuffer) {free (fullBuffer); fullBuffer=NULL; }
		if (otherBuffer) {free (otherBuffer); otherBuffer=NULL; }
		printf ("#Windows %d, #IndelWindows %d maxWindowSize %d numIndels %d\n", curWindow, indelWindowCount, maxWindowSize, numIndels);
		free (posTable);
		posTable = NULL;
		if (indelPositions != NULL)
		{
			sprintf (suffix, "%s_%s.indel-0-Pos", vcf_file, fileList[i].fileName);
			FILE *fpIndels = fopen (suffix, "w");
			assert (fpIndels != NULL);
			fwrite (indelPositions, sizeof(int), numIndels, fpIndels);
			fclose (fpIndels);
			free (indelPositions);
			indelPositions=NULL;
			indelArraySize=0;
			numIndels=0;
		}
		free (ref_chars[i]);
		ref_chars[i]=NULL;
        	struct timeval et;
	        gettimeofday(&et, NULL);
	}*/
	for (i=0; i<NUM_CONTIGS; i++)
	{
		free (ref_chars[i]);
	}
	free (ref_chars);
	free (si);
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
