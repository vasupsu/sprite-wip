#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

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

int main(int argc, char **argv)
{
	fullRec fr;
	otherRec Or;
	int i=0;
	if (argc <3)
	{
		printf ("Usage: decode <SAM file (full or other)> <0-full 1-other>\n");
		return 0;
	}
	int opt = atoi(argv[2]);

	FILE *fp = fopen(argv[1], "r");
	assert(fp!=NULL);

	switch(opt)
	{
		case 0:
			fread (&fr, sizeof(fullRec), 1, fp);
			printf ("#Pos, flag, Qual, MatchLen, seq\n");
			while (!feof(fp))
			{
				printf ("%u, %u, %u, %u, ",  fr.pos, fr.flag, fr.qual, fr.matchLen);
				for (i=0; i<fr.matchLen; i++)
				{
					int curInd=i/4;
					int subInd= i&3;
					uint8_t chInt=(fr.seq[curInd]>>((3-subInd)*2))&3;
					printf ("%c", "ACTG"[chInt]);
				}
				printf ("\n");
				fread (&fr, sizeof(fullRec), 1, fp);
			}
			fclose(fp);
			break;
		case 1:
			fread (&Or, sizeof(otherRec), 1, fp);
			printf ("#Pos, flag, Qual, nCigarOps, seq\t CIGAR\n");
			while (!feof(fp))
			{
				printf ("%u, %u, %u, %u - ",  Or.pos, Or.flag, Or.qual, Or.n_cigar);
				if (Or.n_cigar > 10) Or.n_cigar=10;
				for (i=0; i<Or.n_cigar; i++)
				{
					char c="MIDSH"[Or.cigar[i] & 0xF];
					uint32_t len = Or.cigar[i]>>4;
					printf ("%d%c,", len, c);
				}
				for (i=0; i<101; i++)
				{
					int curInd=i/4;
					int subInd= i&3;
					uint8_t chInt=(Or.seq[curInd]>>((3-subInd)*2))&3;
					printf ("%c", "ACTG"[chInt]);
				}
				printf ("\n");
				fread (&Or, sizeof(otherRec), 1, fp);
			}
			fclose(fp);
			break;
		default:
			printf ("incorrect option\n");
			break;
	}
	return 0;
}
