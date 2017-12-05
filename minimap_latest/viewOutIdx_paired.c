#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>
#include <sys/stat.h>

typedef struct {
        size_t offset;
        size_t offset2;
	int fileNo;
}fileIndex;

struct stat buf;

int  main(int argc, char **argv)
{
	if (argc != 2) {
		fprintf(stderr, "%s fqIdxFile\n", argv[0]);
		exit(1);
	}
	int status = stat (argv[1], &buf);
        assert (status == 0);
	int numRecs = buf.st_size/sizeof(fileIndex);
	int i=0,j=0;
	fileIndex *fI = (fileIndex *)malloc(numRecs * sizeof(fileIndex));
	assert (fI != NULL);
	FILE *fp = fopen (argv[1], "r");
	fread(fI, sizeof(fileIndex), numRecs, fp);
	fclose (fp);
	uint64_t numKmers=0;
	for (i=0; i<numRecs; i++)
	{
		printf ("%d %d %lu,%lu\n", i, fI[i].fileNo, fI[i].offset, fI[i].offset2);
/*        uint64_t chunkKmers = 0;
        for (j=0; j<1048576; j++)
        {
            chunkKmers += fI[i].kmerFreqCount[j];
        }
        printf ("%lu\n", chunkKmers);
        numKmers += chunkKmers;*/
	}
//    printf ("TotalKmers %lu\n", numKmers);
	free (fI);
	return 0;
}
