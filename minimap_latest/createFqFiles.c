#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <sys/stat.h>

int main()
{
	char buf[1000];
	long i=0;

	FILE *fpFq = fopen ("/scratch/03201/tg824998/Venter/aeb_sprite-4/minimap_op_0/unmapped.fq", "r");
	assert (fpFq != NULL);
	fgets (buf, 300, fpFq);
	while (!feof (fpFq))
	{
		int len = strlen(buf);
		buf[len-3]=0;
		printf ("%s\n", buf);
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
