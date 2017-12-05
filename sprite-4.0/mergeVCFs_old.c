#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int main(int argc, char **argv)
{
	FILE *fpParsnip, *fpStrelka;
	if (argc != 4)
	{
		printf ("Usage: a.out <parsnip vcf> <strelka vcf> <output vcf filename>\n");
		return 1;
	}
	char *parsnip_line = (char *)malloc (2000);
	assert (parsnip_line != NULL);
	char *strelka_line = (char *)malloc (200000000);
	assert (strelka_line != NULL);

	char strelka_contigstr[200], posName[100], refBase[100], altBase[100], filter[200];
	int strelka_pos, strelka_qual;
	char parsnip_contigstr[200];
	int parsnip_pos;

	fpParsnip = fopen (argv[1], "r");
	assert (fpParsnip != NULL);
	fgets (parsnip_line, 2000, fpParsnip);
	while (!feof (fpParsnip) && (parsnip_line[0]=='#'))
	{
		fgets (parsnip_line, 2000, fpParsnip);
	}
	if (!feof (fpParsnip))
	{
		sscanf (parsnip_line, "%s\t%d", parsnip_contigstr, &parsnip_pos);
	}
	fpStrelka = fopen (argv[2], "r");
	assert (fpStrelka != NULL);
	FILE *outVcf = fopen (argv[3], "w");
	assert (outVcf != NULL);

	fgets (strelka_line, 200000000, fpStrelka);
	while (!feof (fpStrelka))
	{
		if (strelka_line[0] != '#')
		{
			if (strcmp (filter, "PASS")==0)
			{
				while (!feof (fpParsnip) && (strcmp(parsnip_contigstr, strelka_contigstr) || (parsnip_pos < strelka_pos)))
				{
					fprintf (outVcf, "%s", parsnip_line);
					fgets (parsnip_line, 2000, fpParsnip);
					if (!feof (fpParsnip))
						sscanf (parsnip_line, "%s\t%d", parsnip_contigstr, &parsnip_pos);
				}
				fprintf (outVcf, "%s", strelka_line);
				if ((!feof (fpParsnip)) && !strcmp(parsnip_contigstr, strelka_contigstr) && (parsnip_pos == strelka_pos))
				{
					fgets (parsnip_line, 2000, fpParsnip);
					if (!feof (fpParsnip))
						sscanf (parsnip_line, "%s\t%d", parsnip_contigstr, &parsnip_pos);
				}
			}
		}
		else
			fprintf (outVcf, "%s", strelka_line);
		fgets (strelka_line, 200000000, fpStrelka);
		if (!feof (fpStrelka))
		{
			if (strelka_line[0] != '#')
				sscanf (strelka_line, "%s\t%d\t%s\t%s\t%s\t%d\t%s", strelka_contigstr, &strelka_pos, posName, refBase, altBase, &strelka_qual, filter);
		}
	}
	fclose (fpParsnip);
	fclose (fpStrelka);
	fclose (outVcf);

	free (parsnip_line);
	free (strelka_line);
	return 0;
}
