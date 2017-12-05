#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef USE_OMP
#include <omp.h>
#endif
#include <stdint.h>
#include <sys/stat.h>
#include <assert.h>
#include <sys/time.h>

int bwa_P;
char sam_prefix[500];
int NUM_CONTIGS=93;
int NUM_FILES=0;
int max_out_files = 0;
int numOutChunks, numMaxChunks, maxChunkSize;
char fai_file[300];
int numtasks=1, rank=0, numThreads = 16;
typedef struct 
{
	char *fileName;
	long numRecords;
	long seqLength;
}fileInfo;
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
/*typedef struct
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
*/
fileInfo *fileList=NULL;
long maxSeqLen=0;
int startContig=-1, endContig=-1;
int *procAssign=NULL;

static int sizecmpfn(const void *a, const void *b) {
    fileInfo *av = (fileInfo *) a;
    fileInfo *bv = (fileInfo *) b;
    return av->numRecords-bv->numRecords;
}
void split_fasta (int numChunks, fileInfo *anns, int numContigs, int *numOutChunks, int *numMaxChunks, int *maxChunkSize)
{
        int32_t maxContigSize = 0;
        int i=0;
        for (i=0; i<numContigs; i++)
        {
                if (anns[i].seqLength > maxContigSize) maxContigSize = anns[i].seqLength;
        }
        int curChunks = 0, curContigSize=maxContigSize;
        while (curChunks < numChunks)
        {
                curChunks = 0;
                for (i=0; i<numContigs; i++)
                {
                        curChunks += (anns[i].seqLength/curContigSize)+(anns[i].seqLength%curContigSize > 0);
                }
		if (rank == 0)
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
                int curChunksPerContig = (anns[i].seqLength/curContigSize)+(anns[i].seqLength%curContigSize > 0);
                curChunks += curChunksPerContig;
                if (curChunksPerContig > maxChunksPerContig) maxChunksPerContig = curChunksPerContig;
        }
        *numOutChunks = curChunks;
        *numMaxChunks = maxChunksPerContig;
        *maxChunkSize = curContigSize;
}
int main(int argc, char **argv)
{
	struct stat buf;
	int i=0,j=0,k=0;
	if (argc != 5)
	{
		printf ("Args: <maxOutFiles> <ref.fasta.fai> <ain,aeb file prefix> <#BWA MPI processes>\n");
                printf ("maxOutFiles - Max number of output files per process\nref.fasta.fai - Absolute path of the fai file corresponding to the input fasta file.\naeb,aib file prefix - prefix of the filename of aeb and aib files created by Prune.\n#BWA MPI processes - Number of MPI tasks for PRUNE run\n");
		return 0;
	}
	max_out_files = atoi (argv[1]);
	strcpy(fai_file, argv[2]);
	FILE *fp_fai = fopen(fai_file, "r");
	assert(fp_fai != NULL);
	strcpy (sam_prefix, argv[3]);
	bwa_P = atoi(argv[4]);
	
#ifdef USE_MPI
	MPI_Init(&argc,&argv);
        MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
	if (rank == 0)
		printf ("NumTasks %d rank %d\n", numtasks, rank);
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
		contigNo++;
		fgets(line, 200, fp_fai);
        }
	fclose(fp_fai);
	NUM_CONTIGS = contigNo;
	split_fasta (max_out_files, fileList, NUM_CONTIGS, &numOutChunks, &numMaxChunks, &maxChunkSize);
	if (rank == 0)
		printf ("NUM_CONTIGS %d numOutChunks %d numMaxChunks %d maxChunkSize %d\n", NUM_CONTIGS, numOutChunks, numMaxChunks, maxChunkSize);
	size_t *nrecs = (size_t *) calloc (numMaxChunks*NUM_CONTIGS, sizeof(size_t));
	assert (nrecs != NULL);
	for (i=0; i<NUM_CONTIGS; i++)
	{
		for (j=0; j<numMaxChunks; j++)
		{
			for (k=0; k<bwa_P; k++)
			{
				sprintf(fileName, "%s_%d/C%d_%d.aeb", sam_prefix, k, i, j); 
				int status=stat (fileName, &buf);
				if (status != -1)
				{
//					printf ("%s\n", fileName);
					nrecs[i*numMaxChunks+j] += (buf.st_size/sizeof(fullRec));
				}

				sprintf(fileName, "%s_%d/unmapped_%d_%d.aeb", sam_prefix, k, i, j); 
				status=stat (fileName, &buf);
				if (status != -1)
				{
//					printf ("%s\n", fileName);
					nrecs[i*numMaxChunks+j] += (buf.st_size/sizeof(fullRec));
				}

				sprintf(fileName, "%s_%d/unmapped_%d_%d.aib", sam_prefix, k, i, j); 
				status=stat (fileName, &buf);
				if (status != -1)
				{
//					printf ("%s\n", fileName);
					nrecs[i*numMaxChunks+j] += (buf.st_size/sizeof(otherRec));
				}
			}
			if ((rank == 0) && (nrecs[i*numMaxChunks+j] > 0))
				printf ("[%d] Contig %d Segment %d - %d recs\n", i*numMaxChunks+j, i, j, nrecs[i*numMaxChunks+j]);
			totalRecords += nrecs[i*numMaxChunks+j];
		}
	}
	if (rank == 0)
		printf ("TotalRecords %ld\n", totalRecords);
//	qsort(fileList, NUM_CONTIGS, sizeof(fileInfo), sizecmpfn);
	long recordsPerProc = totalRecords / numtasks;
	procAssign = (int *)malloc(sizeof(int) * (NUM_CONTIGS*numMaxChunks + 1));
	bzero(procAssign, sizeof(int) * (numtasks+1));
	procAssign[numtasks]=NUM_CONTIGS*numMaxChunks;
	long prefixSum = nrecs[0];
	int index=0,fileInd=0;
	for (index=1; index<numtasks; index++)
	{
		while (prefixSum < index*recordsPerProc)
		{
			fileInd++;
			prefixSum += nrecs[fileInd];
		}
		procAssign[index] = fileInd;
	}
	if (rank == 0)
	{
		for (i=0; i <= numtasks; i++)
		{
			printf ("%d\t", procAssign[i]);
		}
		printf ("\n");
	}


	
//	for (i=procAssign[rank]; i<procAssign[rank+1]; )
	for (contigNo = 0; contigNo < NUM_CONTIGS; contigNo++)
	{
			////////////////count #records aligning to each position////////////
		int segNo;
#ifdef USE_OMP
		#pragma omp parallel for num_threads(numThreads) private(segNo, i, j, k, fileName) schedule(dynamic)
#endif
		for (segNo=0; segNo<numMaxChunks; segNo++)
		{

#ifdef USE_OMP
			int tid = omp_get_thread_num();
#else
			int tid = 0;
#endif
			struct timeval stime, etime;
			gettimeofday (&stime, NULL);
			fullRec *Fr=NULL;
			otherRec *Or=NULL;
	
			i = contigNo * numMaxChunks + segNo;
//			if (i!=33) continue;
			if ((i<procAssign[rank]) || (i>=procAssign[rank+1])) continue;
			int segSize = (maxChunkSize < fileList[contigNo].seqLength) ? maxChunkSize : fileList[contigNo].seqLength;
			int segStart = segNo*segSize;
			if ((segStart+segSize) >= fileList[contigNo].seqLength) segSize = fileList[contigNo].seqLength - segStart;
			if (segSize <= 0) continue;

			int *contigPosCount=NULL;
			contigPosCount = (int *)malloc(segSize * sizeof(int));
			assert (contigPosCount != NULL);
			bzero (contigPosCount, segSize * sizeof(int));
			int numAeb=0, numAib=0;
			int totalAlnRecs = 0;
//			printf ("Rank %d tid %d Contig(%d) Segment %d Length %d, startPos %d\n", rank, tid, contigNo, segNo, segSize, segStart);
			Fr = (fullRec *)malloc(100000*sizeof(fullRec));
	        	Or = (otherRec *)malloc(100000*sizeof(otherRec));	
			for (j=0; j<bwa_P; j++)
			{
				sprintf (fileName, "%s_%d/C%d_%d.aeb", sam_prefix, j, contigNo, segNo);
				FILE *fp_sam = fopen(fileName, "r");
	                        if(fp_sam != NULL)
				{
					int nread=fread (Fr, sizeof(fullRec), 100000, fp_sam);
					while (!feof(fp_sam))
					{
						for (k=0; k<nread; k++)
						{
		//					printf ("1:k %d pos %d\n", k, Fr[k].pos-1);
							if ((Fr[k].pos >= fileList[contigNo].seqLength) || (Fr[k].pos <= 0))
							{
								printf ("2:k %d pos %u\n", k, Fr[k].pos-1);
								continue;
							}
							contigPosCount[Fr[k].pos-1-segStart]++;
							totalAlnRecs++;
						}
						nread=fread (Fr, sizeof(fullRec), 100000, fp_sam);
					}
				
					for (k=0; k<nread; k++)
					{
		//				printf ("3:k %d pos %u\n", k, Fr[k].pos-1);
						if ((Fr[k].pos >= fileList[contigNo].seqLength) || (Fr[k].pos <= 0))
						{
							printf ("4:k %d pos %u\n", k, Fr[k].pos-1);
							continue;
						}
						contigPosCount[Fr[k].pos-1-segStart]++;
						totalAlnRecs++;
					}
					fclose(fp_sam);
				}
				sprintf (fileName, "%s_%d/unmapped_%d_%d.aeb", sam_prefix, j, contigNo, segNo);
				fp_sam = fopen(fileName, "r");
	                        if(fp_sam != NULL)
				{
					int nread=fread (Fr, sizeof(fullRec), 100000, fp_sam);
					while (!feof(fp_sam))
					{
						for (k=0; k<nread; k++)
						{
		//					printf ("1:k %d pos %d\n", k, Fr[k].pos-1);
							if ((Fr[k].pos >= fileList[contigNo].seqLength) || (Fr[k].pos <= 0))
							{
								printf ("2:k %d pos %u\n", k, Fr[k].pos-1);
								continue;
							}
							contigPosCount[Fr[k].pos-1-segStart]++;
							totalAlnRecs++;
						}
						nread=fread (Fr, sizeof(fullRec), 100000, fp_sam);
					}
					
					for (k=0; k<nread; k++)
					{
		//				printf ("3:k %d pos %u\n", k, Fr[k].pos-1);
						if ((Fr[k].pos >= fileList[contigNo].seqLength) || (Fr[k].pos <= 0))
						{
							printf ("4:k %d pos %u\n", k, Fr[k].pos-1);
							continue;
						}
						contigPosCount[Fr[k].pos-1-segStart]++;
						totalAlnRecs++;
					}
					fclose(fp_sam);
				}
			}
			int *startPos=NULL;
			startPos = (int *)malloc(segSize * sizeof(int));
			assert (startPos != NULL);
			startPos[0]=0;
			for (j=1; j<segSize; j++)
			{
				startPos[j]=startPos[j-1]+contigPosCount[j-1];
			}
			numAeb = totalAlnRecs;
//			printf ("[%d,%d] numAEBrecs %d\n", contigNo, segNo, startPos[segSize-1]);
			fullRec *frWrite = (fullRec *)malloc(totalAlnRecs * sizeof(fullRec));
			assert(frWrite != NULL);
			for (j=0; j<bwa_P; j++)
			{
				sprintf (fileName, "%s_%d/C%d_%d.aeb", sam_prefix, j, contigNo, segNo);
				FILE *fp_sam = fopen(fileName, "r");
	                        if(fp_sam != NULL)
				{
					int nread=fread (Fr, sizeof(fullRec), 100000, fp_sam);
					while (!feof(fp_sam))
					{
						for (k=0; k<nread; k++)
						{
		//					printf ("5:k %d pos %u\n", k, Fr[k].pos-1);
							if ((Fr[k].pos >= fileList[contigNo].seqLength) || (Fr[k].pos <= 0))
							{
								printf ("6:k %d pos %u\n", k, Fr[k].pos-1);
								continue;
							}
							memcpy(&(frWrite[startPos[Fr[k].pos-1-segStart]++]), &(Fr[k]), sizeof(fullRec));
						}
						nread=fread (Fr, sizeof(fullRec), 100000, fp_sam);
					}
					
					for (k=0; k<nread; k++)
					{
		//				printf ("7:k %d pos %u\n", k, Fr[k].pos-1);
						if ((Fr[k].pos >= fileList[contigNo].seqLength) || (Fr[k].pos <= 0))
						{
							printf ("8:k %d pos %u\n", k, Fr[k].pos-1);
							continue;
						}
						memcpy(&(frWrite[startPos[Fr[k].pos-1-segStart]++]), &(Fr[k]), sizeof(fullRec));
					}
					fclose(fp_sam);
		//			int ret=unlink(fileName);
		//                        if (ret==-1)
		//                                printf ("File %s Not Deleted. Error (%d)\n", fileName, errno);
				}
				sprintf (fileName, "%s_%d/unmapped_%d_%d.aeb", sam_prefix, j, contigNo, segNo);
				fp_sam = fopen(fileName, "r");
	                        if(fp_sam != NULL)
				{
					int nread=fread (Fr, sizeof(fullRec), 100000, fp_sam);
					while (!feof(fp_sam))
					{
						for (k=0; k<nread; k++)
						{
		//					printf ("5:k %d pos %u\n", k, Fr[k].pos-1);
							if ((Fr[k].pos >= fileList[contigNo].seqLength) || (Fr[k].pos <= 0))
							{
								printf ("6:k %d pos %u\n", k, Fr[k].pos-1);
								continue;
							}
							memcpy(&(frWrite[startPos[Fr[k].pos-1-segStart]++]), &(Fr[k]), sizeof(fullRec));
						}
						nread=fread (Fr, sizeof(fullRec), 100000, fp_sam);
					}
					
					for (k=0; k<nread; k++)
					{
		//				printf ("7:k %d pos %u\n", k, Fr[k].pos-1);
						if ((Fr[k].pos >= fileList[contigNo].seqLength) || (Fr[k].pos <= 0))
						{
							printf ("8:k %d pos %u\n", k, Fr[k].pos-1);
							continue;
						}
						memcpy(&(frWrite[startPos[Fr[k].pos-1-segStart]++]), &(Fr[k]), sizeof(fullRec));
					}
					fclose(fp_sam);
		//			int ret=unlink(fileName);
		//                        if (ret==-1)
		//                                printf ("File %s Not Deleted. Error (%d)\n", fileName, errno);
				}
			}
	
			if (totalAlnRecs > 0)
			{
				sprintf (fileName, "%s_0/C%d_%d_sorted.aeb", sam_prefix, contigNo, segNo);
				FILE *fp_out = fopen(fileName, "w");
				assert(fp_out!=NULL);
				fwrite (frWrite, sizeof(fullRec), totalAlnRecs, fp_out);
				fclose (fp_out);
			}
	////////////////////////////////end of full match
			free(frWrite);
			totalAlnRecs = 0;
			bzero (contigPosCount, segSize * sizeof(int));
			bzero (startPos, segSize * sizeof(int));
			////////////count other recs
			for (j=0; j<bwa_P; j++)
			{
				sprintf (fileName, "%s_%d/unmapped_%d_%d.aib", sam_prefix, j, contigNo, segNo);
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
						contigPosCount[Or[k].pos-1-segStart]++;
						totalAlnRecs++;
					}
					nread=fread (Or, sizeof(otherRec), 100000, fp_sam);
				}
				
				for (k=0; k<nread; k++)
				{
					contigPosCount[Or[k].pos-1-segStart]++;
					totalAlnRecs++;
				}
				fclose(fp_sam);
			}
			startPos[0]=0;
			//File Offset for each ref position
	                for (j=1; j<segSize; j++)
	                {
	                        startPos[j]=startPos[j-1]+contigPosCount[j-1];
	                }
//			printf ("[%d,%d]numAIBrecs %d\n", contigNo, segNo, startPos[segSize-1]);
			numAib = totalAlnRecs;
			otherRec *orWrite = (otherRec *)malloc(totalAlnRecs * sizeof(otherRec));
	                assert(orWrite != NULL);
			for (j=0; j<bwa_P; j++)
	                {
	                        sprintf (fileName, "%s_%d/unmapped_%d_%d.aib", sam_prefix, j, contigNo, segNo);
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
	                                        memcpy(&(orWrite[startPos[Or[k].pos-1-segStart]++]), &(Or[k]), sizeof(otherRec));
	                                }
	                                nread=fread (Or, sizeof(otherRec), 100000, fp_sam);
	                        }
	
	                        for (k=0; k<nread; k++)
	                        {
	                                memcpy(&(orWrite[startPos[Or[k].pos-1-segStart]++]), &(Or[k]), sizeof(otherRec));
	                        }
	                        fclose(fp_sam);
	//			int ret=unlink(fileName);
	//                        if (ret==-1)
	//                                printf ("File %s Not Deleted\n", fileName);
	                }
			if (totalAlnRecs > 0)
			{
				sprintf (fileName, "%s_0/C%d_%d_sorted.aib", sam_prefix, contigNo, segNo);
		                FILE *fp_out = fopen(fileName, "w");
	        	        assert(fp_out!=NULL);
	                	fwrite (orWrite, sizeof(otherRec), totalAlnRecs, fp_out);
		                fclose (fp_out);
			}
			free(orWrite);
			
			free (contigPosCount);
			free (startPos);
			free(Fr);
			free(Or);
	
			gettimeofday (&etime, NULL);
			long elapsed = (etime.tv_sec * 1000000 + etime.tv_usec) - (stime.tv_sec * 1000000 + stime.tv_usec);
			if ((numAeb+numAib) > 0)
				printf ("%d %d %d - %lf sec\n", contigNo, segNo, numAeb+numAib, (double)elapsed/1000000);
		}//end of segment
	}//end of contig
	free (nrecs);
#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	for (i=0; i<NUM_CONTIGS; i++)
		free(fileList[i].fileName);
	free(fileList);
	free(procAssign);
#ifdef USE_MPI
	MPI_Finalize();
#endif
	return 0;
}
