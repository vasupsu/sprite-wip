#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

int main(int argc, char **argv)
{
	if (argc != 3)
	{
		printf ("Usage: a.out VCFfile BEDfile\n");
		return 1;
	}
	FILE *bedIn = fopen (argv[2], "r");
	assert (bedIn != NULL);
	FILE *vcfIn = fopen (argv[1], "r");
	assert (vcfIn != NULL);
	char bed_line[500], vcf_line[50000];
	fgets (bed_line, 500, bedIn);
	char chrName[50];
	int sPos, ePos;
	sscanf (bed_line, "%s\t%d\t%d", chrName, &sPos, &ePos);
	fgets (vcf_line, 50000, vcfIn);
	while (!feof (vcfIn) && vcf_line[0]=='#')
	{
		printf ("%s", vcf_line);
		fgets (vcf_line, 50000, vcfIn);
	}
	char vcfChrName[50];
	int vcfPos;
	sscanf (vcf_line, "%s\t%d", vcfChrName, &vcfPos);
	while (!feof (bedIn))
	{
		int vcfChanged=0;
//		fprintf (stdout, "B%s:%s", bed_line, vcf_line);
		while (!feof (vcfIn) && (strcmp (chrName, vcfChrName) || (vcfPos < sPos)))
		{
			char oldVcfChrName[50];
			strcpy (oldVcfChrName, vcfChrName);
			fgets (vcf_line, 50000, vcfIn);
			if (!feof (vcfIn))
			{
				sscanf (vcf_line, "%s\t%d", vcfChrName, &vcfPos);
				if (!strcmp (oldVcfChrName, chrName) && strcmp(oldVcfChrName, vcfChrName))
				{
					vcfChanged=1;
					break;
				}
			}
		}
//		fprintf (stdout, "A%s:%s\n", bed_line, vcf_line);
		while (!feof (vcfIn) && !strcmp (chrName, vcfChrName) && (vcfPos < ePos))
		{
			printf ("%s", vcf_line);
			fgets (vcf_line, 50000, vcfIn);
			char oldVcfChrName[50];
			strcpy (oldVcfChrName, vcfChrName);
			if (!feof (vcfIn))
			{
				sscanf (vcf_line, "%s\t%d", vcfChrName, &vcfPos);
				if (strcmp(oldVcfChrName, vcfChrName))
					vcfChanged=1;
			}
		}
		fgets (bed_line, 500, bedIn);
		if (!feof (bedIn))
		{
			sscanf (bed_line, "%s\t%d\t%d", chrName, &sPos, &ePos);
			if (vcfChanged)
			{
//				fprintf (stderr, "%s:%d-%d vcfchanged %d vcfChrName %s\n", chrName, sPos, ePos, vcfChanged, vcfChrName);
				while (!feof (bedIn) && strcmp(vcfChrName, chrName))
				{
					fgets (bed_line, 500, bedIn);
					if (!feof (bedIn))
						sscanf (bed_line, "%s\t%d\t%d", chrName, &sPos, &ePos);
				}
			}
		}
	}

	fclose (bedIn);
	fclose (vcfIn);
	return 0;
}
