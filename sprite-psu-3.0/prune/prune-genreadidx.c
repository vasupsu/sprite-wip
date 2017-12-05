#include <stdio.h>
#include <zlib.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>

int chunkSize=100000;//reads
int main(int argc, char **argv)
{
	int numtasks, rank=0;
	numtasks=1;
	gzFile gz_fp;
	int gzFd;
	if (argc != 3)
	{
		printf ("Usage: a.out read.fq read1Length\n");
		return 0;
	}
	struct stat buf;
	int i=0,j;
	
	char *file1 = argv[1];

	int status = stat (argv[1], &buf);
	assert(status == 0);

        long fileLen = buf.st_size;
 
	gzFd = open(file1, O_RDONLY);
        assert (gzFd >= 0);
        gz_fp = gzdopen(gzFd, "r");
	
	int read1Len=atoi(argv[2]);
	char *readChunk = (char *)malloc(read1Len*8);

	long writeOfs = 0;
	FILE *fpOp = NULL;

	long startOfs=(rank*chunkSize*read1Len*4);
	if (rank==0)
		startOfs=numtasks*chunkSize*read1Len*4;
	long endOfs =(chunkSize*read1Len*4);
	if (rank == (numtasks-1))
		endOfs=fileLen;
	int count=1;
	char opFile[200];
	
	int eof=0;
	sprintf (opFile, "%s.idx", argv[1]);
	printf ("IndexFile:%s\n", opFile);
	fpOp=fopen(opFile, "w");
	if (rank==0)
	fwrite (&writeOfs, sizeof(long), 1, fpOp);
	assert (fpOp != NULL);
	while (eof==0)
	{
		gzseek(gz_fp, (z_off_t)startOfs, SEEK_SET);
		size_t nread =gzread (gz_fp, readChunk, read1Len* 8);
		if (nread != (read1Len* 8))
			eof=1;
		if (eof) break;
		int flag=0;
		for (i=0; i<nread; i++)
		{
			if ((readChunk[i]=='\n') && (readChunk[i+1]=='+') && (readChunk[i+2]=='\n'))
			{
				flag=1;
				i=i+3;
			}
			if (flag && (readChunk[i]=='\n'))
			{
				i++;
//				for (j=i; readChunk[j]!='\n'; j++)
//					printf ("%c", readChunk[j]);
//				printf ("\n");
//				break;
				for (j=i; readChunk[j]!='/'; j++);
				if (readChunk[j+1]=='2')
				{
					flag=0;
					i=j;
				}
				else
					break;
			}
		}
		assert((flag==1) || (eof==1));
		writeOfs = startOfs+i;
		fwrite (&writeOfs, sizeof(long), 1, fpOp);
		printf ("writeOfs %ld\n", writeOfs);
		count ++;
		startOfs+=numtasks*chunkSize*read1Len*4;
	}
	printf ("Created %d chunks\n", count);
	free (readChunk);
	gzclose (gz_fp);
	if (fpOp != NULL)
		fclose (fpOp);
	return 0;
}
