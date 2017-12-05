#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>

int main(int argc, char **argv)
{
	char line[1000];
	FILE *fp = fopen (argv[1], "r");
	assert (fp != NULL);
	fseek (fp, atol (argv[2]), SEEK_SET);
	fread (line, 1, 998, fp);
	line[999]='\0';
	printf ("%s\n", line);
	fclose (fp);
	return 0;
}
