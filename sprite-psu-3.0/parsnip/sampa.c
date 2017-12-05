#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#ifdef USE_MPI
#include <mpi.h>
#endif
#include <stdint.h>
#include <sys/stat.h>
#include <assert.h>

int bwa_P;
char sam_prefix[500];
int NUM_CONTIGS=93;
int NUM_FILES=0;
char fai_file[300];
typedef struct 
{
	char *fileName;
	long numRecords;
	long seqLength;
}fileInfo;
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

fileInfo *fileList=NULL;
long maxSeqLen=0;
int *recCount=NULL;
int startContig=-1, endContig=-1;
int *procAssign=NULL;
fullRec *Fr=NULL;
otherRec *Or=NULL;
int *contigPosCount=NULL;
int *startPos=NULL;

static int sizecmpfn(const void *a, const void *b) {
    fileInfo *av = (fileInfo *) a;
    fileInfo *bv = (fileInfo *) b;
    return av->numRecords-bv->numRecords;
}
int main(int argc, char **argv)
{
	struct stat buf;
	int i=0,j=0,k=0;
	if (argc < 4)
	{
		printf ("Args: <path-to-fasta-index .fai> <ain,aeb file prefix> <#BWA MPI processes>\n");
		printf ("Args: <ref.fasta.fai> <aeb,aib file prefix> <#BWA MPI processes>\n");
                printf ("ref.fasta.fai - Absolute path of the fai file corresponding to the input fasta file.\naeb,aib file prefix - prefix of the filename of aeb and aib files created by Prune.\n#BWA MPI processes - Number of MPI tasks for PRUNE run\n");
		return 0;
	}
	strcpy(fai_file, argv[1]);
	FILE *fp_fai = fopen(fai_file, "r");
	assert(fp_fai != NULL);
	strcpy (sam_prefix, argv[2]);
	bwa_P = atoi(argv[3]);
	int numtasks=1, rank=0;
	
#ifdef USE_MPI
	MPI_Init(&argc,&argv);
        MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif

	char line[200];
	NUM_CONTIGS=0;
        fgets(line, 200, fp_fai);
        while (!feof(fp_fai))
        {
		NUM_CONTIGS++;
		fgets(line, 200, fp_fai);
        }
	fclose(fp_fai);

	fp_fai = fopen(fai_file, "r");
        long length=0;
	fileList = (fileInfo *)malloc(NUM_CONTIGS * sizeof(fileInfo));
        fgets(line, 200, fp_fai);
        int contigNo=0;
	char contigName[100];
	char fileName[500];
	long totalRecords=0;
        while (!feof(fp_fai))
        {
		sscanf(line, "%s\t%ld", contigName, &length);
		int fLen = strlen(contigName);
		fileList[contigNo].fileName = (char *)malloc(fLen+1);
		strcpy(fileList[contigNo].fileName, contigName);
		fileList[contigNo].numRecords = 0;
		fileList[contigNo].seqLength = length;
		if (maxSeqLen < length)
			maxSeqLen = length;
		for (i=0; i<bwa_P; i++)
		{
			sprintf(fileName, "%s_%s_%d.aeb", sam_prefix, contigName, i); 
			int status=stat (fileName, &buf);
			if (status != -1)
				fileList[contigNo].numRecords += (buf.st_size/sizeof(fullRec));
		}
		for (i=0; i<bwa_P; i++)
		{
			sprintf(fileName, "%s_%s_%d.aib", sam_prefix, contigName, i); 
			int status=stat (fileName, &buf);
			if (status != -1)
				fileList[contigNo].numRecords += (buf.st_size/sizeof(otherRec));
		}
		totalRecords += fileList[contigNo].numRecords;
		contigNo++;
		fgets(line, 200, fp_fai);
        }
	fclose(fp_fai);
	NUM_CONTIGS = contigNo;
//	qsort(fileList, NUM_CONTIGS, sizeof(fileInfo), sizecmpfn);
	recCount = (int *)malloc(sizeof(int) * maxSeqLen);
	assert(recCount != NULL);
	long recordsPerProc = totalRecords / numtasks;
	procAssign = (int *)malloc(sizeof(int) * (numtasks+1));
	bzero(procAssign, sizeof(int) * (numtasks+1));
	procAssign[numtasks]=NUM_CONTIGS;
	long prefixSum=fileList[0].numRecords;
	int index=0,fileInd=0;
	for (index=1; index<numtasks; index++)
	{
		while (prefixSum < index*recordsPerProc)
		{
			fileInd++;
			prefixSum += fileList[fileInd].numRecords;
		}
		procAssign[index] = fileInd;
	}

	Fr = (fullRec *)malloc(100000*sizeof(fullRec));
        Or = (otherRec *)malloc(100000*sizeof(otherRec));	

	
	for (i=procAssign[rank]; i<procAssign[rank+1]; i++)
//	for (i=22; i<23; i++)
	{
		printf ("Rank %d Contig(%d) %s\n", rank, i, fileList[i].fileName);
		int totalAlnRecs = 0;
		contigPosCount = (int *)malloc(fileList[i].seqLength * sizeof(int));
		assert (contigPosCount != NULL);
		bzero (contigPosCount, fileList[i].seqLength * sizeof(int));
		////////////////count #records aligning to each position////////////
		for (j=0; j<bwa_P; j++)
		{
			sprintf (fileName, "%s_%s_%d.aeb", sam_prefix, fileList[i].fileName, j);
			FILE *fp_sam = fopen(fileName, "r");
                        if(fp_sam == NULL)
			{
				continue;
			}
			int nread=fread (Fr, sizeof(fullRec), 100000, fp_sam);
			while (!feof(fp_sam))
			{
				for (k=0; k<nread; k++)
				{
					contigPosCount[Fr[k].pos-1]++;
					totalAlnRecs++;
				}
				nread=fread (Fr, sizeof(fullRec), 100000, fp_sam);
			}
			
			for (k=0; k<nread; k++)
			{
				contigPosCount[Fr[k].pos-1]++;
				totalAlnRecs++;
			}
			fclose(fp_sam);
		}
		startPos = (int *)malloc(fileList[i].seqLength * sizeof(int));
		assert (startPos != NULL);
		startPos[0]=0;
		for (j=1; j<fileList[i].seqLength; j++)
		{
			startPos[j]=startPos[j-1]+contigPosCount[j-1];
		}
		fullRec *frWrite = (fullRec *)malloc(totalAlnRecs * sizeof(fullRec));
		assert(frWrite != NULL);
		for (j=0; j<bwa_P; j++)
		{
			sprintf (fileName, "%s_%s_%d.aeb", sam_prefix, fileList[i].fileName, j);
			FILE *fp_sam = fopen(fileName, "r");
                        if(fp_sam == NULL)
			{
				continue;
			}
			int nread=fread (Fr, sizeof(fullRec), 100000, fp_sam);
			while (!feof(fp_sam))
			{
				for (k=0; k<nread; k++)
				{
					memcpy(&(frWrite[startPos[Fr[k].pos-1]++]), &(Fr[k]), sizeof(fullRec));
				}
				nread=fread (Fr, sizeof(fullRec), 100000, fp_sam);
			}
			
			for (k=0; k<nread; k++)
			{
				memcpy(&(frWrite[startPos[Fr[k].pos-1]++]), &(Fr[k]), sizeof(fullRec));
			}
			fclose(fp_sam);
//			int ret=unlink(fileName);
//                      if (ret==-1)
//                                printf ("File %s Not Deleted. Error (%d)\n", fileName, errno);
		}

		if (totalAlnRecs > 0)
		{
			sprintf (fileName, "%s_%s_sorted.aeb", sam_prefix, fileList[i].fileName);
			FILE *fp_out = fopen(fileName, "w");
			assert(fp_out!=NULL);
			fwrite (frWrite, sizeof(fullRec), totalAlnRecs, fp_out);
			fclose (fp_out);
		}
////////////////////////////////end of full match
		free(frWrite);
		totalAlnRecs = 0;
		bzero (contigPosCount, fileList[i].seqLength * sizeof(int));
		bzero (startPos, fileList[i].seqLength * sizeof(int));
		////////////count other recs
		for (j=0; j<bwa_P; j++)
		{
			sprintf (fileName, "%s_%s_%d.aib", sam_prefix, fileList[i].fileName, j);
			FILE *fp_sam = fopen(fileName, "r");
                        if(fp_sam == NULL)
			{
				continue;
			}
			int nread=fread (Or, sizeof(otherRec), 100000, fp_sam);
			while (!feof(fp_sam))
			{
				for (k=0; k<nread; k++)
				{
					contigPosCount[Or[k].pos-1]++;
					totalAlnRecs++;
				}
				nread=fread (Or, sizeof(otherRec), 100000, fp_sam);
			}
			
			for (k=0; k<nread; k++)
			{
				contigPosCount[Or[k].pos-1]++;
				totalAlnRecs++;
			}
			fclose(fp_sam);
		}
		startPos[0]=0;
		//File Offset for each ref position
                for (j=1; j<fileList[i].seqLength; j++)
                {
                        startPos[j]=startPos[j-1]+contigPosCount[j-1];
                }
		otherRec *orWrite = (otherRec *)malloc(totalAlnRecs * sizeof(otherRec));
                assert(orWrite != NULL);
		for (j=0; j<bwa_P; j++)
                {
                        sprintf (fileName, "%s_%s_%d.aib", sam_prefix, fileList[i].fileName, j);
                        FILE *fp_sam = fopen(fileName, "r");
                        if(fp_sam == NULL)
			{
				continue;
			}
                        int nread=fread (Or, sizeof(otherRec), 100000, fp_sam);
			while (!feof(fp_sam))
                        {
                                for (k=0; k<nread; k++)
                                {
                                        memcpy(&(orWrite[startPos[Or[k].pos-1]++]), &(Or[k]), sizeof(otherRec));
                                }
                                nread=fread (Or, sizeof(otherRec), 100000, fp_sam);
                        }

                        for (k=0; k<nread; k++)
                        {
                                memcpy(&(orWrite[startPos[Or[k].pos-1]++]), &(Or[k]), sizeof(otherRec));
                        }
                        fclose(fp_sam);
//			int ret=unlink(fileName);
//                        if (ret==-1)
//                                printf ("File %s Not Deleted\n", fileName);
                }
		if (totalAlnRecs > 0)
		{
			sprintf (fileName, "%s_%s_sorted.aib", sam_prefix, fileList[i].fileName);
	                FILE *fp_out = fopen(fileName, "w");
        	        assert(fp_out!=NULL);
                	fwrite (orWrite, sizeof(otherRec), totalAlnRecs, fp_out);
	                fclose (fp_out);
		}
		free(orWrite);
		
		free (contigPosCount);
		free (startPos);
	}
#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	for (i=0; i<NUM_CONTIGS; i++)
		free(fileList[i].fileName);
	free(fileList);
	free(recCount);
	free(procAssign);
	free(Fr);
	free(Or);
#ifdef USE_MPI
	MPI_Finalize();
#endif
	return 0;
}
