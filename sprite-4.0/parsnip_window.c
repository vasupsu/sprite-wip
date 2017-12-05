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
int WINDOW_SIZE=50000;
int maxBufferRecords=1000000;
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

char *ref_chars=NULL;
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
ref_pos *posTable=NULL;

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
        uint32_t pos;
	uint16_t flag;
        uint8_t seq[60];
	uint8_t quals[120];
        uint8_t qual;
        uint8_t n_cigar;
        uint16_t cigar[10];
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

int *indelPositions=NULL;
int numIndels=0;
int indelArraySize=0;

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
            ref_chars[ind]=ch;
        }
	free (contigStr);
}
inline int processOther (otherRec *rec, int winStart)
{	
	int ret=0;
	int i=0,j=0;
	int pos=0;
	int startPos = rec->pos-1;
	startPos = startPos % WINDOW_SIZE;
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
						if (ref_chars[startPos+winStart] != Seq[j])
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
				if (numIndels==indelArraySize)
				{
					indelArraySize=(indelArraySize==0)?1000:indelArraySize<<1;
					indelPositions = (int *)realloc (indelPositions, indelArraySize*sizeof(int));
					assert (indelPositions != NULL);
				}
				indelPositions[numIndels++]=startPos+winStart;
				pos+=len;
				break;
			case 2://Delete
				ret=1;
				if (numIndels==indelArraySize)
				{
					indelArraySize=(indelArraySize==0)?1000:indelArraySize<<1;
					indelPositions = (int *)realloc (indelPositions, indelArraySize*sizeof(int));
					assert (indelPositions != NULL);
				}
				indelPositions[numIndels++]=startPos+winStart;
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
	return ret;
}
void outputVcfRecs(int startPos, FILE *fp_vcf, ref_pos *pos_table, char *cName, int seqLen, FILE *fn){
	long ind=0;
//	fprintf (fp_vcf, "Window: %d-%d\n", startPos, startPos + seqLen-1);
	for (ind=0; ind<WINDOW_SIZE; ind++)
        {
        	char ch=ref_chars[ind+startPos];
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
                if (((startPos+ind+1) == 891059) || ((startPos+ind+1) == 1900232) || ((startPos+ind+1) == 1335422))
                {
			int avg_mapq = 0;
			if (max_per > ref_per)
				avg_mapq = (pos_table[ind].count[4]==0)?0:(pos_table[ind].mapQ/(pos_table[ind].count[4]-ref_per));
	                printf("%s\t%ld\t.\t%c\t%c\t%d\t.", cName, startPos+ind+1, ch, max_char, avg_mapq);
        	        printf("\tDP=%d;ACTG=%d,%d,%d,%d;FC=%u;RC=%u;FT=%u;RT=%u;AAF=%.2f\n", pos_table[ind].count[4], pos_table[ind].count[0], pos_table[ind].count[1], pos_table[ind].count[2], pos_table[ind].count[3], pos_table[ind].fwdCount, pos_table[ind].revCount, pos_table[ind].fwdTotal, pos_table[ind].revTotal, maxPer);
                }
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
		                   fprintf(fp_vcf, "%s\t%ld\t.\t%c\t%c\t%d\t.", cName, startPos+ind+1, ch, max_char, avg_mapq);
		                   fprintf(fp_vcf, "\tDP=%d;ACTG=%d,%d,%d,%d;FC=%u;RC=%u;FT=%u;RT=%u;SB=%.2f;AAF=%.2f", pos_table[ind].count[4], pos_table[ind].count[0], pos_table[ind].count[1], pos_table[ind].count[2], pos_table[ind].count[3], pos_table[ind].fwdCount, pos_table[ind].revCount, pos_table[ind].fwdTotal, pos_table[ind].revTotal, r1, maxPer);
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
        readFasta (contig, fp_fasta, posTable, seqLen);
	assert ((WINDOW_SIZE%10)==0);
}
void parsnip_aeb(fullRec *Fr, int winStart)
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
		pos = pos % WINDOW_SIZE;
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
				if (Seq[l] != ref_chars[pos+winStart])
        	                {
                	                posTable[pos].mapQ+=Fr->qual;
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
int parsnip_aib(otherRec *Or, int winStart)
{
        totalAlnFiles=1;
	return processOther(Or, winStart);
}
void parsnip_write(char *cName, int seqLen, int startPos, FILE *fp_vcf)
{
        if (totalAlnFiles > 0)
        {
                assert(fp_vcf!=NULL);
                outputVcfRecs(startPos, fp_vcf, posTable, cName, seqLen, NULL);
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
//	for (i=/*procAssign[rank]*/; i<procAssign[rank+1]; i++)
	for (i=1; i<2; i++)
	{
		printf ("Contig %d\n", i);
        	struct timeval tm;
	        gettimeofday(&tm, NULL);

		posTable = (ref_pos *)malloc((WINDOW_SIZE+WINDOW_SIZE/10) * sizeof(ref_pos));
		ref_chars = (char *) malloc (fileList[i].seqLength * sizeof (char));
		assert (ref_chars != NULL);
		printf ("Alloc postable for contig %d - %ld entries\n", i, WINDOW_SIZE+WINDOW_SIZE/10/*fileList[i].seqLength*/);
		assert(posTable!=NULL);
		bzero(posTable, (WINDOW_SIZE+WINDOW_SIZE/10)/*fileList[i].seqLength*/ * sizeof(ref_pos));

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
			parsnip_init(fileList[i].seqLength,fasta_file, i);
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
/*			if (!initDone || (Or->pos < winStart))
			{
				winStart = Or->pos-1;
				winEnd = Or->pos + seqLen - 2;
			}
			if ((Or->pos < winEnd) && ((Or->pos + seqLen - 2) > winEnd))
				winEnd = Or->pos + seqLen - 2;*/
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
        	char suffix[500];
                sprintf (suffix, "%s.vcf", vcf_file);
                FILE *fp_vcf = fopen(suffix, "w");
 		while ((fp_aeb && (curAebRec < maxAebRecs)) || (fp_aib && (curAibRec < maxAibRecs)))
		{
			if (fp_aeb && (curAebRec < maxAebRecs))
			{
				Fr = &fullBuffer[curAebRec];
				if ((fp_aib && (fp_aib && (curAibRec < maxAibRecs)) && (Fr->pos <= otherBuffer[curAibRec].pos)) || (!fp_aib) || (fp_aib && (curAibRec >= maxAibRecs)))
				{
					while (Fr->pos > winEnd)
					{
//						if (!indelWindow)
						{//writeWindow
			                        	parsnip_write(fileList[i].fileName, winEnd-winStart, winStart, fp_vcf);
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
						parsnip_aeb (Fr, winStart);
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
                        				parsnip_write(fileList[i].fileName, winEnd-winStart, winStart, fp_vcf);
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
						indelWindow |= parsnip_aib (Or, winStart);
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
			else
			{

				Or = &otherBuffer[curAibRec];
				int seqLen = getSeqLen(Or);
				while (Or->pos > winEnd)
				{
//					if (!indelWindow)
					{//writeWindow
						parsnip_write(fileList[i].fileName, winEnd-winStart, winStart, fp_vcf);
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
					indelWindow |= parsnip_aib (Or, winStart);
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
//                parsnip_write(fileList[i].fileName, fileList[i].seqLength/*winEnd-winStart+1*/, 0/*winStart*/, fp_vcf);
                parsnip_write(fileList[i].fileName, winEnd-winStart, winStart, fp_vcf);
		bzero (posTable, WINDOW_SIZE*sizeof(ref_pos));
		memmove (posTable, &posTable[WINDOW_SIZE], (WINDOW_SIZE/10)*sizeof(ref_pos));
		bzero (&posTable[WINDOW_SIZE], (WINDOW_SIZE/10)*sizeof(ref_pos));
		winStart = winEnd;
		winEnd = winStart + WINDOW_SIZE;
                parsnip_write(fileList[i].fileName, winEnd-winStart, winStart, fp_vcf);
		if (fp_aeb) fclose (fp_aeb);
		if (fp_aib) fclose (fp_aib);
/*		free (Fr);
		Fr = NULL;
		free (Or);
		Or = NULL;*/
		if (fullBuffer) {free (fullBuffer); fullBuffer=NULL; }
		if (otherBuffer) {free (otherBuffer); otherBuffer=NULL; }
		fclose (fp_vcf);
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
		free (ref_chars);
		ref_chars=NULL;
        	struct timeval et;
	        gettimeofday(&et, NULL);
	}
#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
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
	MPI_Finalize();
#endif
	return 0;
}
