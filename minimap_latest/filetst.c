#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>

int main()
{
	char buf[10000];
	FILE *fp = fopen ("in.paf", "r");
	assert (fp != NULL);
	size_t nread = fread (buf, 1, 10000, fp);
	printf ("nread %lu buf[%lu]=%c(%d)\n", nread, nread-1, buf[nread-1], buf[nread-1]);
	fclose (fp);
	return 0;
}
