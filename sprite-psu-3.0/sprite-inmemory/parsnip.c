#include <stdio.h>
#include <sys/time.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#ifdef USE_MPI
#include <mpi.h>
#endif
#include <math.h>
#include <stdint.h>
#include <sys/stat.h>
#include <assert.h>

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
typedef struct ref_pos
{
   uint32_t count[5];   
   uint32_t readPos[4];
   uint16_t fwdCount;//#alt alleles in fwd reads
   uint16_t fwdTotal;//#total fwd reads
   uint16_t revCount;
   uint16_t revTotal;
   char ref_char;
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
ref_pos *posTable=NULL;

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

long maxSeqLen=0;
int startContig=-1, endContig=-1;
int *procAssign=NULL;
int *startPos=NULL;
int *contigPosCount=NULL;

typedef struct parsnip_data
{
        fullRec *fR;
        long len;
        int contig;
        int tid;
        long *nextOfs;
        long *remLen;
        long *startTime;
        long *endTime;
}parsnip_data;

void readFasta (int contigNo, FILE *fp_fasta, ref_pos *posTable, int len, long startPos)
{
	long ind=0;
	fseek(fp_fasta, startPos, SEEK_SET);
	int lenToRead = len + (len/50)+100;
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
            posTable[ind].ref_char=ch;
        }
	free (contigStr);
}
inline void processOther (otherRec *rec, ref_pos *posTable)
{	
	int i=0,j=0;
	int pos=0;
	int startPos = rec->pos-1;
	int ind=0;
	char Seq[110];
	if (rec->n_cigar > 10) rec->n_cigar=10;
	for (i=0; i<rec->n_cigar; i++)
	{
		int cigarType=rec->cigar[i] & 0xF;
		uint32_t len = rec->cigar[i]>>4;
		switch(cigarType)
		{
			case 0://Match
				for(j=0;j<len;j++)
				{
					posTable[startPos].count[4]++;
					int curInd=pos/4;
				        int subInd= pos&3;
                               		ind=(rec->seq[curInd]>>((3-subInd)*2))&3;
					Seq[j]="ACTG"[ind];
						
					posTable[startPos].count[ind]++;
					if (rec->flag & 0x10)
						posTable[startPos].revTotal++;
					else
						posTable[startPos].fwdTotal++;
					posTable[startPos].readPos[ind]+=pos;
					if (posTable[startPos].ref_char != Seq[j])
					{
						posTable[startPos].mapQ+=rec->qual;
						if (rec->flag & 0x10)
							posTable[startPos].revCount++;
						else
							posTable[startPos].fwdCount++;
					}
					pos++;
					startPos++;
				}
				Seq[j]=0;
				break;
			case 1://Insert
				pos+=len;
				break;
			case 2://Delete
				for(j=0;j<len;j++)
				{
					posTable[startPos].count[4]++;
					startPos++;
				}
				break;
			case 3://Soft
				pos+=len;
				break;
			case 4://H
//				pos+=len;
				break;
		}
	}
}
void outputVcfRecs(int startPos, FILE *fp_vcf, ref_pos *pos_table, char *cName, int seqLen, FILE *fn){
	long ind=0;
	for (ind=startPos; ind<(startPos+seqLen); ind++)
        {
        	char ch=pos_table[ind].ref_char;
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
                if ((maxPer > 20) && (max_char != ch) && ((max_per > 1)) && (pos_table[ind].count[4] > 0))
                {
		   if (maxPer > 80) gt1=1;
//				   float fwdRatio = (pos_table[ind].fwdTotal>0)?((float)(pos_table[ind].fwdCount)/(pos_table[ind].fwdTotal)):0;
//				   float revRatio = (pos_table[ind].revTotal>0)?((float)(pos_table[ind].revCount)/(pos_table[ind].revTotal)):0;
				   int isZero = 0;
				   if ((pos_table[ind].fwdTotal ==0) || (pos_table[ind].revTotal ==0)) isZero=1;
//				   float ratioOfRatio = (fwdRatio<revRatio)?(fwdRatio/revRatio):(revRatio/fwdRatio);
				   float r1 = (pos_table[ind].fwdCount < pos_table[ind].revCount)?((float)pos_table[ind].fwdCount/pos_table[ind].revCount):((float)pos_table[ind].revCount/pos_table[ind].fwdCount);
				int avg_mapq = (pos_table[ind].count[4]==0)?0:(pos_table[ind].mapQ/(pos_table[ind].count[4]-ref_per));
				if ( ((r1>.1) || ((r1==0) && (isZero==1))) && (avg_mapq >= 20))
				{	
		                   fprintf(fp_vcf, "%s\t%ld\t.\t%c\t%c\t%d\t.", cName, ind+1, ch, max_char, avg_mapq);
		                   fprintf(fp_vcf, "\tDP=%d;ACTG=%d,%d,%d,%d;FC=%u;RC=%u;FT=%u;RT=%u;SB=%.2f;AAF=", pos_table[ind].count[4], pos_table[ind].count[0], pos_table[ind].count[1], pos_table[ind].count[2], pos_table[ind].count[3], pos_table[ind].fwdCount, pos_table[ind].revCount, pos_table[ind].fwdTotal, pos_table[ind].revTotal, r1);
				   if (maxPer < 25) fprintf(fp_vcf, "20");
				   else if (maxPer <30) fprintf(fp_vcf, "25");
				   else if (maxPer <35) fprintf(fp_vcf, "30");
				   else fprintf(fp_vcf, "35");
		                   fprintf(fp_vcf, "\tGT\t%d|1\n", gt1);
				}
                }
        }
}
void parsnip_init(int rank, int seqLen, char *fasta_file, char *vcf_prefix, int numContigs, int contig, int nt)
{
        char line[200];
        long length=0,start=0;
	char fai_file[200];
	char contigName[200];
	sprintf (fai_file, "%s.fai", fasta_file);
        FILE *fp_fai = fopen(fai_file, "r");
	assert (fp_fai!=NULL);
        fgets(line, 200, fp_fai);
        int contigNo=0;
        while (!feof(fp_fai) && (contigNo <= contig))
        {
                sscanf(line, "%s\t%ld\t%ld", contigName, &length, &start);
                contigNo++;
                fgets(line, 200, fp_fai);
        }
        fclose(fp_fai);
	assert (contigNo == (contig+1));

        posTable = (ref_pos *)malloc(seqLen * sizeof(ref_pos));
        assert(posTable != NULL);
        bzero(posTable, seqLen * sizeof(ref_pos));
        sprintf (fileName, "%s_c%d", vcf_prefix, contig);
        fp_fasta = fopen(fasta_file, "r");
	assert (fp_fasta != NULL);
        readFasta (contig, fp_fasta, posTable, seqLen, start);
}
void *parsnip_aeb(void *data)
{
        parsnip_data *pData = (parsnip_data *)data;
        fullRec *Fr = pData->fR;
        long numFr = pData->len;
        int j, l;
        char Seq[150];

        long numRecords1 = numFr;
        if (numRecords1 > 0)totalAlnFiles++;
        for (j=0; j<numRecords1; j++)
        {
                unsigned int pos=Fr[j].pos-1;
                if ((pData->nextOfs!=NULL) && (pData->tid!=(numThreads-1)) &&((pos + 200) >= Fr[*(pData->nextOfs)].pos))
                {
                        *(pData->nextOfs)=j;
                        *(pData->remLen)=numFr-j;
                        break;
                }
                for (l=0; l<Fr[j].matchLen; l++)
                {
                        int curInd=l/4;
                        int subInd= l&3;
                        uint8_t ind=(Fr[j].seq[curInd]>>((3-subInd)*2))&3;
			Seq[l]="ACTG"[ind];
                        posTable[pos].count[ind]++;
                        posTable[pos].count[4]++;
			if (Fr[j].flag & 0x10)
                                posTable[pos].revTotal++;
                        else
                                posTable[pos].fwdTotal++;
			if (Seq[l] != posTable[pos].ref_char)
                        {
                                posTable[pos].mapQ+=Fr[j].qual;
                                if (Fr[j].flag & 0x10)
                                     	posTable[pos].revCount++;
                                else
                                        posTable[pos].fwdCount++;
                        }
			pos++;
		}
		Seq[l]=0;
	}
        return NULL;
}			
void parsnip_aib(otherRec *Or, long numOr, int contig)
{
        int j;
        long numRecords2 = numOr;
        if (numRecords2 > 0)totalAlnFiles++;
        for (j=0; j<numRecords2; j++)
        {
                       processOther(&Or[j], posTable);
        }
}
void *parsnip_write(void *data)
{
        parsnip_data *pData = (parsnip_data *)data;
        char *cName = (char *)pData->fR;
        int seqLen = (int)pData->len;
        int startPos = pData->contig;
        int tid=pData->tid;
        char suffix[500];
        if (totalAlnFiles > 0)
        {
                sprintf (suffix, "%s_t%d.vcf", fileName, tid);
                FILE *fp_vcf = fopen(suffix, "w");
                assert(fp_vcf!=NULL);
                outputVcfRecs(startPos, fp_vcf, posTable, cName, seqLen, NULL);
                fclose(fp_vcf);
        }
	return NULL;
}
