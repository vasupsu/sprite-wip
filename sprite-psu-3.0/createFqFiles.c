#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/stat.h>
int main()
{
	char buf[1000];
	long i=0;

	FILE *fpFq = fopen ("/scratch/03201/tg824998/Venter/reads1.fq", "r");
	assert (fpFq != NULL);
	fgets (buf, 300, fpFq);
	while (!feof (fpFq))
	{
		printf ("@Venter.%ld/1\n", i);
		fgets (buf, 300, fpFq);
		printf ("%s", buf);
		fgets (buf, 300, fpFq);
		printf ("%s", buf);
		fgets (buf, 300, fpFq);
		printf ("%s", buf);
		fgets (buf, 300, fpFq);
		i++;
//		if (i==10)
//			break;
	}
	fclose (fpFq);
	return 0;
}
