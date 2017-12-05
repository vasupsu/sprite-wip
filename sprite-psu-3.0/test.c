#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/stat.h>
int main()
{
	struct stat buf1;
	int ret = stat("/scratch/03201/tg824998/Venter/read1.fq.idx", &buf1);
	assert (ret == 0);
	int numChunks = buf1.st_size/sizeof(long);
	printf ("numChunks %d\n", numChunks);
	long *idx = (long *)malloc (numChunks *sizeof(long));
	assert (idx != NULL);

	FILE *fpIdx = fopen ("/scratch/03201/tg824998/Venter/read1.fq.idx", "r");
	assert (fpIdx != NULL);
	fread (idx, sizeof(long), numChunks, fpIdx);
	fclose (fpIdx);

	char buf[1000];
	int i=0;

	FILE *fpFq = fopen ("/scratch/03201/tg824998/Venter/read2.fq", "r");
	assert (fpFq != NULL);
	for (i=0; i<numChunks; i++)
	{
		long ofs = idx[i];
		fseek (fpFq, ofs, SEEK_SET);
		fread (buf, 1, 1000, fpFq);
		int j=0;
		printf ("[%d] Ofs %ld: ", i, ofs);
		while (buf[j] != '\n')
		{
			printf ("%c", buf[j]);
			j++;
			if (j>200) break;
		}
		printf ("\n");
	}
	fclose (fpFq);
	free (idx);
	return 0;
}
