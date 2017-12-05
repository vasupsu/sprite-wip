#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>

#define CHUNK_SIZE 200000000

uint32_t getFirstRead (char *buf)
{
	int i=0;
	int found=0;
	for (;i<1000;i++)
	{
		if ((buf[i]=='\n') && (buf[i+1]=='+') && (buf[i+2]=='\n'))
		{
			i+=3;
			while (buf[i]!='\n') i++;
			while (buf[i] != '.') i++;
			i++;
			int i1=i;
			while (buf[i1] != ' ') 
			{
				assert ((buf[i1] >= '0') && (buf[i1] <= '9'));
				i1++;
			}
			buf[i1]='\0';
//			printf ("%s\n", &buf[i]);
			return (atoi (&buf[i]));
		}
		
	}
	return 0;
}
size_t printRead (char *buf, uint32_t targetReadId)
{
	int i=0;
	int found=0;
	size_t readOfs=0;
	buf[100]='\0';
//	printf ("buf %s\n", buf);
	for (;i<CHUNK_SIZE;i++)
	{
		if ((buf[i]=='\n') && (buf[i+1]=='+') && (buf[i+2]=='\n'))
		{
			i+=3;
			while (buf[i]!='\n') i++;i++;
			readOfs = i;
			int i1=i;
//	buf [i1+100]=0;
			uint32_t curReadId = 0;
			while (buf[i1] != '.') 
			{
//				printf ("%c", buf[i1]);
				i1++;
			}
			i1++;
//	printf ("buf %s\n", &buf[i1]);
			while ((buf[i1] != ' ') && (buf[i1] != '\n'))
			{
				assert ((buf[i1] >= '0') && (buf[i1] <= '9'));
				curReadId = curReadId*10 + buf[i1] - '0';
				i1++;
			}
//			printf ("curReadId %u\n", curReadId);
			if (curReadId == targetReadId)
			{
				found=1;
				break;
			}
		}
	}
	int numLines=0;
	assert (found==1);
	while (numLines < 4)
	{
		if (found)
			printf ("%c", buf[i]);
		if (buf[i]=='\n') numLines++;
		i++;
	}
	return readOfs;
}
static int cmpReadId(const void *p1, const void *p2)
{
        uint32_t f1=*((uint32_t *)p1);
        uint32_t f2=*((uint32_t *)p2);
           return (f1-f2);
}
#define N_TARGETS 33
int main(int argc, char **argv)
{
	if (argc != 2)
	{
		printf ("a.out a.fq\n");
		return 0;
	}
	FILE *fp = fopen (argv[1], "r");
	assert (fp != NULL);

	char buf[1000];
	size_t prevOfs=0, curOfs=CHUNK_SIZE;
	uint32_t targets[N_TARGETS] = {429350674, 429350692, 783258963, 429350702, 429350791, 429350746, 429350761, 429350719, 429350722, 429350693, 429350694, 429350720, 429350695, 429350688, 429350696, 429350717, 429350736, 429350697, 429350698, 429350701, 429350699, 429350700, 429350725, 429350732, 429350865, 429350672, 429350743, 429350726, 429350734, 429350723, 429350724, 429350796, 787265110};
	qsort (targets, N_TARGETS, sizeof(uint32_t), cmpReadId);
	uint32_t targetReadId = targets[0];//1-based
//	printf ("Target %u\n", targetReadId);
	int targetIdx=1;
	size_t startOfs=0;
	fseek (fp, curOfs, SEEK_SET);
	while (!feof (fp))
	{
		size_t nread  = fread (buf, 1, 1000, fp);
//		printf ("nread %d \n", nread);
		if (nread < 1000) break;
		uint32_t firstReadId = getFirstRead (buf);
//		printf ("%u %u\n", firstReadId, targetReadId);
		while (firstReadId > targetReadId)
		{
//			printf ("fseek ofs %lu\n", prevOfs);
			fseek (fp, prevOfs, SEEK_SET);
			char *chunk = (char *) malloc (CHUNK_SIZE);
			assert (chunk != NULL);
			size_t nread = fread (chunk, 1, CHUNK_SIZE, fp);
//			printf ("after fread: nread %lu\n", nread);
			startOfs = printRead (chunk, targetReadId)+prevOfs;
			free (chunk);
			targetReadId = targets[targetIdx++];
		}
		if (targetIdx == N_TARGETS) break;
		prevOfs = curOfs;
		curOfs += CHUNK_SIZE;
		fseek (fp, curOfs, SEEK_SET);
	}
	fclose (fp);

	return 0;
}
