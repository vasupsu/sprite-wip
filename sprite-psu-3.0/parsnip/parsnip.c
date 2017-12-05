#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <errno.h>
#include <ctype.h>
#include <unistd.h>
#ifdef USE_MPI
#include <mpi.h>
#endif
#include <math.h>
#include <getopt.h>
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

static int sizecmpfn(const void *a, const void *b) {
    fileInfo *av = (fileInfo *) a;
    fileInfo *bv = (fileInfo *) b;
    return av->numRecords-bv->numRecords;
}
void readFasta (int contigNo, FILE *fp_fasta, ref_pos *posTable, int len)
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
					if (rec->qual > 0)
						posTable[startPos].count[4]++;
					int curInd=pos/4;
				        int subInd= pos&3;
                               		ind=(rec->seq[curInd]>>((3-subInd)*2))&3;
					Seq[j]="ACTG"[ind];
					if (rec->qual >0)
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
					if (rec->qual >0)
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
//	printf ("outputVcfRecs begin\n");
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
                if ((maxPer > AAF) && (max_char != ch) && ((max_per > AAC)) && (pos_table[ind].count[4] > DP))
                {
		   if (maxPer > 80) gt1=1;
//				   float fwdRatio = (pos_table[ind].fwdTotal>0)?((float)(pos_table[ind].fwdCount)/(pos_table[ind].fwdTotal)):0;
//				   float revRatio = (pos_table[ind].revTotal>0)?((float)(pos_table[ind].revCount)/(pos_table[ind].revTotal)):0;
				   int isZero = 0;
				   if ((pos_table[ind].fwdTotal ==0) || (pos_table[ind].revTotal ==0)) isZero=1;
//				   float ratioOfRatio = (fwdRatio<revRatio)?(fwdRatio/revRatio):(revRatio/fwdRatio);
				   float r1 = (pos_table[ind].fwdCount < pos_table[ind].revCount)?((float)pos_table[ind].fwdCount/pos_table[ind].revCount):((float)pos_table[ind].revCount/pos_table[ind].fwdCount);
				int avg_mapq = (pos_table[ind].count[4]==0)?0:(pos_table[ind].mapQ/(pos_table[ind].count[4]-ref_per));
				if ( ((r1>SB) || ((r1==0) && (isZero==1))) && (avg_mapq >= MQ))
				{	
		                   fprintf(fp_vcf, "%s\t%ld\t.\t%c\t%c\t%d\t.", cName, ind+1, ch, max_char, avg_mapq);
		                   fprintf(fp_vcf, "\tDP=%d;ACTG=%d,%d,%d,%d;FC=%u;RC=%u;FT=%u;RT=%u;SB=%.2f;AAF=", pos_table[ind].count[4], pos_table[ind].count[0], pos_table[ind].count[1], pos_table[ind].count[2], pos_table[ind].count[3], pos_table[ind].fwdCount, pos_table[ind].revCount, pos_table[ind].fwdTotal, pos_table[ind].revTotal, r1);
				   if (maxPer < 25) fprintf(fp_vcf, "20");
				   else if (maxPer <30) fprintf(fp_vcf, "25");
				   else if (maxPer <35) fprintf(fp_vcf, "30");
				   else fprintf(fp_vcf, "35");
		                   fprintf(fp_vcf, "\tGT\t%d/1\n", gt1);
				}
/*				
			}*/
                }
        }
//	printf ("outputVcfRecs end\n");
}
void parsnip_init(int rank, int seqLen, char *fasta_file, char *vcf_prefix, int numContigs, int contig, int nt)
{
/*        posTable = (ref_pos *)malloc(seqLen * sizeof(ref_pos));
        assert(posTable != NULL);
        bzero(posTable, seqLen * sizeof(ref_pos));*/
        sprintf (fileName, "%s_c%d", vcf_file, fileList[contig].contigNo);
        fp_fasta = fopen(fasta_file, "r");
	assert (fp_fasta != NULL);
        readFasta (contig, fp_fasta, posTable, seqLen);
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
                long pos=Fr[j].pos-1;
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
			if (Fr[j].qual > 0)
			{
	                        posTable[pos].count[ind]++;
        	                posTable[pos].count[4]++;
			}
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
                fprintf (stderr, "Usage: parsnip [options] <path to fasta file> <sorted aeb,aib prefix> <output vcf file prefix>\n");
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
	printf ("rank %d/%d numThreads %d\n", rank, numtasks, numThreads);
        char line[500];
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
                sprintf(fileName, "%s_%s_sorted.aeb", sam_prefix, contigName);
                int r1 = stat (fileName, &buf);
                if (r1==0)
		{
                        fileList[contigNo].numRecords += (buf.st_size/sizeof(fullRec));
			fileList[contigNo].numAeb = (buf.st_size/sizeof(fullRec));
		}

                sprintf(fileName, "%s_%s_sorted.aib", sam_prefix, contigName);
                int r2 = stat (fileName, &buf);
                if (r2==0)
		{
                        fileList[contigNo].numRecords += (buf.st_size/sizeof(otherRec));
			fileList[contigNo].numAib = (buf.st_size/sizeof(otherRec));
		}
                totalRecords += fileList[contigNo].numRecords;
                contigNo++;
                fgets(line, 200, fp_fai);
        }
        fclose(fp_fai);
//        qsort(fileList, NUM_CONTIGS, sizeof(fileInfo), sizecmpfn);
        long recordsPerProc = totalRecords / numtasks;
        procAssign = (int *)malloc(sizeof(int) * (numtasks+1));
        bzero(procAssign, sizeof(int) * (numtasks+1));
        procAssign[numtasks]=NUM_CONTIGS;
        long prefixSum=fileList[0].numRecords;
        int fileInd=0;
        long index=0;

	for (index=1; index<numtasks; index++)
        {
                while (prefixSum < index*recordsPerProc)
                {
                        fileInd++;
                        prefixSum += fileList[fileInd].numRecords;
                }
                procAssign[index] = fileInd;
        }
	printf ("Rank %d: Contigs [%d,%d)\n", rank, procAssign[rank], procAssign[rank+1]);
	for (i=procAssign[rank]; i<26; i++)
//	for (i=22; i<26; i++)
	{
		printf ("[%d]Contig %d\n", rank, i);
        	struct timeval tm;
	        gettimeofday(&tm, NULL);

		posTable = (ref_pos *)calloc(fileList[i].seqLength, sizeof(ref_pos));
		assert(posTable!=NULL);
//		printf ("Alloc postable for contig %d - %ld entries (%ld bytes)\n", i, fileList[i].seqLength, fileList[i].seqLength*sizeof(ref_pos));
//		bzero(posTable, fileList[i].seqLength * sizeof(ref_pos));

		int initDone=0;
		totalAlnFiles=0;
		sprintf(fileName, "%s_%s_sorted.aeb", sam_prefix, fileList[i].fileName);
		parsnip_data *pData=(parsnip_data *)malloc(numThreads*sizeof(parsnip_data));
                pthread_t *thr= (pthread_t *)malloc(numThreads*sizeof(pthread_t));
       	        assert(pData!=NULL);
		FILE *fp_sam = fopen(fileName, "r");
		if (fp_sam != NULL)
		{
			fullRec *Fr = (fullRec *)malloc(fileList[i].numAeb * sizeof(fullRec));
	        	assert(Fr != NULL);
			long nRead = fread (Fr, sizeof(fullRec), fileList[i].numAeb, fp_sam);
			assert (nRead == fileList[i].numAeb);
			fclose (fp_sam);
//			printf ("before parsnip init\n");
			parsnip_init(rank,fileList[i].seqLength,fasta_file, "parsnip", NUM_CONTIGS, i, numThreads);
			initDone = 1;
                	long *nextOfs=(long *)malloc(numThreads*sizeof(long));
	                long *remLen=(long *)malloc(numThreads*sizeof(long));
        	        long *aebStart=(long *)malloc(numThreads*sizeof(long));
                	assert((nextOfs != NULL) && (remLen != NULL) && (aebStart !=NULL));
	                for (j=0; j<numThreads; j++)
        	        {
                	        pthread_attr_t attr;
                        	pthread_attr_init(&attr);
	                        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
				aebStart[j]=fileList[i].numAeb/numThreads*j;
				long len = fileList[i].numAeb/numThreads*(j+1);
				len = (j==(numThreads-1))?fileList[i].numAeb:len;
				len -= aebStart[j];
				nextOfs[j]=len;
        	                remLen[j]=0;
                	        pData[j].fR=&(Fr[aebStart[j]]);
				pData[j].len=len;
	                        pData[j].contig=i;
        	                pData[j].tid=j;
                	        pData[j].nextOfs = &(nextOfs[j]);
                        	pData[j].remLen = &(remLen[j]);
//				printf ("create therad %d\n", j);
	                        pthread_create(&thr[j], &attr, parsnip_aeb, &pData[j]);
			}
			for (j=0; j<numThreads; j++)
	                {
        	                pthread_join(thr[j], 0);
//				printf ("join thread %d\n", j);
                        	pData[j].fR=&(Fr[aebStart[j]+nextOfs[j]]);
	                        pData[j].len=remLen[j];
        	                pData[j].nextOfs=NULL;
                	        pData[j].remLen=NULL;
				if (pData[j].len > 0)
	                        	parsnip_aeb(&pData[j]);
//				printf ("thread %d after AEB\n", j);
	                }
			free (nextOfs);
			free (remLen);
			free (aebStart);
			free (Fr);
			Fr = NULL;
		}
		else
		{
			printf ("File %s not found\n", fileName);
		}
		sprintf(fileName, "%s_%s_sorted.aib", sam_prefix, fileList[i].fileName);
		fp_sam = fopen(fileName, "r");
		if (fp_sam != NULL)
                {
        		otherRec *Or = (otherRec *)malloc(fileList[i].numAib * sizeof(otherRec));
	        	assert(Or != NULL);
			long nRead = fread (Or, sizeof(otherRec), fileList[i].numAib, fp_sam);
			assert (nRead == fileList[i].numAib);
			fclose (fp_sam);
			if (initDone==0)
			{
				parsnip_init(rank,fileList[i].seqLength,fasta_file, "parsnip", NUM_CONTIGS, i, numThreads);
				initDone=1;
			}
			parsnip_aib(Or, fileList[i].numAib, i);
			free (Or);
			Or = NULL;
		}
		else
		{
			printf ("File %s not found\n", fileName);
		}
		if (initDone)
		{
			for (j=0; j<numThreads; j++)
                	{
				pthread_attr_t attr;
	                        pthread_attr_init(&attr);
        	                pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
                	        long seqStart=fileList[i].seqLength/numThreads*j;
                        	long len = fileList[i].seqLength/numThreads*(j+1);

	                        len = (j==(numThreads-1))?fileList[i].seqLength:len;
        	                len -= seqStart;
	                        pData[j].fR=(fullRec *)fileList[i].fileName;
	                        pData[j].len=len;
        	                pData[j].contig=(int)seqStart;
                	        pData[j].tid=j;
                        	pthread_create(&thr[j], &attr, parsnip_write, &pData[j]);
			}
			for (j=0; j<numThreads; j++)
                        {
				pthread_join(thr[j], 0);
			}
		}
		free (pData);
		free (posTable);
		posTable = NULL;
		free (thr);
        	struct timeval et;
	        gettimeofday(&et, NULL);
	}
#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	char suffix[500];
	char vcfBuf[10000000];
	if (rank==0)
	{
		sprintf (suffix, "%s.vcf", vcf_file);
		FILE *fp_vcf = fopen(suffix, "w");
		assert (fp_vcf != NULL);
                fprintf (fp_vcf, "##fileformat=VCFv4.1\n##source=SPRITE3.0\n##reference=%s\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth at the locus\">\n##INFO=<ID=ACTG,Number=4,Type=Integer,Description=\"Count for each allele type\">\n##INFO=<ID=FC,Number=1,Type=Integer,Description=\"Alternate allele count in forward strands\">\n##INFO=<ID=RC,Number=1,Type=Integer,Description=\"Alternate allele count in reverse strands\">\n##INFO=<ID=FT,Number=1,Type=Integer,Description=\"Total forward strands\">\n##INFO=<ID=RT,Number=1,Type=Integer,Description=\"Total reverse strands\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tvcf\n",fasta_file);
		int nread=0;
		for (i=0; i<NUM_CONTIGS; i++)
		{
			for (j=0; j<numThreads; j++)
			{
				sprintf (suffix, "%s_c%d_t%d.vcf", vcf_file, i, j);
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
					int stat=unlink(suffix);
					if (stat==-1)
						printf ("Cannot remove temporary VCF file %s. Error code %d\n", suffix, errno);
				}
			}
		}
		fclose (fp_vcf);
	}
#ifdef USE_MPI
	MPI_Finalize();
#endif
	return 0;
}
