#include "sam.h"
#include "hts.h"
#include "bgzf.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <omp.h>

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

typedef struct
{
   FILE *fp_fullMatch;
   FILE *fp_others;
   char *fullFileName;
   char *otherFileName;
   long curFullCount;
   long curOtherCount;
   long outputChunkSize;
   fullRec *Fr;
   otherRec *Or;
}contigInfo;

contigInfo *auxInfo=NULL;

void readBamRec (bam1_t *b, bam_hdr_t *h, hts_itr_t *iter, htsFile *fp, FILE *fp_aeb, FILE *fp_aib)
{
	int i=0;
	bam1_core_t *c;
	char *data;
	fullRec FR;
        otherRec OR;

//	int id=omp_get_thread_num();
	int curTid=-1;
	long totUnaligned=0;
	int ret=0;
	int count=0;
	while (ret >= 0)
	{
		c = &b->core;
		data = (char *)b->data;
		int *cigar=(int *)(data+c->l_qname);
		uint8_t *seq = bam_get_seq(b);
//		uint8_t *qualStr = bam_get_qual(b);
		if (c->flag & 0x4) 
		{
			totUnaligned++;
			ret=sam_itr_next(fp, iter, b);
			continue;
		}
		if ((curTid==-1) || (curTid!=c->tid))
		{
			if (curTid!=-1) break;
			curTid = c->tid;
		}
		
		if ((c->n_cigar==1) && (cigar[0]%16==0))
                {
                        char ch;
                        uint8_t chInt;
			if (auxInfo[i].fp_fullMatch == NULL)
		                auxInfo[i].fp_fullMatch = fopen(auxInfo[i].fullFileName, "w");
                        FR.pos = c->pos+1;
                        FR.qual = c->qual;
			FR.flag = c->flag;
                        FR.matchLen = cigar[0]/16;
                        for (i = 0; i < c->l_qseq; ++i)
                        {
                                ch="=ACMGRSVTWYHKDBN"[bam_seqi(seq, i)];
                                chInt = (uint8_t)ch & 0xf;
				int curInd=i/4;
				switch (i&3)
                                {       
                                        case 0: 
                                                FR.seq[curInd]=(ch<<5)&0xc0;
                                                break;
                                        case 1: 
                                                FR.seq[curInd] |= ((ch<<3)&0x30);
                                                break;
                                        case 2: 
                                                FR.seq[curInd] |= ((ch<<1)&0xc);
                                                break;
                                        case 3: 
                                                FR.seq[curInd] |= ((ch>>1) & 0x3);
                                                break;
                                }
                        }
                        if (auxInfo[c->tid].curFullCount == auxInfo[c->tid].outputChunkSize)
                        {
                                if (auxInfo[c->tid].fp_fullMatch == NULL)
                                {
                                        auxInfo[c->tid].fp_fullMatch = fopen(auxInfo[c->tid].fullFileName, "w");
                                        assert(auxInfo[c->tid].fp_fullMatch != NULL);
                                }
                                fwrite(auxInfo[c->tid].Fr, sizeof(fullRec), auxInfo[c->tid].curFullCount, auxInfo[c->tid].fp_fullMatch);
                                auxInfo[c->tid].curFullCount = 0;
                        }
                        memcpy(&(auxInfo[c->tid].Fr[auxInfo[c->tid].curFullCount]), &FR, sizeof(fullRec));
                        auxInfo[c->tid].curFullCount++;
                }
		else
                {
			if (auxInfo[i].fp_others == NULL)
                		auxInfo[i].fp_others = fopen(auxInfo[i].otherFileName, "w");
                        char ch;
                        uint8_t chInt;
                        OR.pos = c->pos+1;
			OR.flag = c->flag;
                        OR.qual = c->qual;
                        OR.n_cigar = c->n_cigar;
			for (i=0; i<c->n_cigar; i++)
			{
				OR.cigar[i]=(uint16_t)cigar[i];
				if (i > 9) break;
			}

                        for (i = 0; i < c->l_qseq; ++i)
                        {
                                ch="=ACMGRSVTWYHKDBN"[bam_seqi(seq, i)];
                                chInt = (uint8_t)ch & 0xf;
				int curInd=i/4;
				switch (i&3)
                                {       
                                        case 0: 
                                                OR.seq[curInd]=(ch<<5)&0xc0;
                                                break;
                                        case 1: 
                                                OR.seq[curInd] |= ((ch<<3)&0x30);
                                                break;
                                        case 2: 
                                                OR.seq[curInd] |= ((ch<<1)&0xc);
                                                break;
                                        case 3: 
                                                OR.seq[curInd] |= ((ch>>1) & 0x3);
                                                break;
                                }
                        }
                        if (auxInfo[c->tid].curOtherCount == auxInfo[c->tid].outputChunkSize)
                        {
                                if (auxInfo[c->tid].fp_others == NULL)
                                {
                                        auxInfo[c->tid].fp_others = fopen(auxInfo[c->tid].otherFileName, "w");
                                        assert(auxInfo[c->tid].fp_others != NULL);
                                }
                                fwrite(auxInfo[c->tid].Or, sizeof(otherRec), auxInfo[c->tid].curOtherCount, auxInfo[c->tid].fp_others);
                                auxInfo[c->tid].curOtherCount = 0;
                        }
                        memcpy(&(auxInfo[c->tid].Or[auxInfo[c->tid].curOtherCount]), &OR, sizeof(otherRec));
                        auxInfo[c->tid].curOtherCount++;
                }
		count ++;
		ret = sam_itr_next(fp, iter, b);
	}
	if (auxInfo[c->tid].curFullCount > 0)
        {
                if (auxInfo[c->tid].fp_fullMatch == NULL)
                {
	                auxInfo[c->tid].fp_fullMatch = fopen(auxInfo[c->tid].fullFileName, "w");
                        assert(auxInfo[c->tid].fp_fullMatch != NULL);
                }
                fwrite(auxInfo[c->tid].Fr, sizeof(fullRec), auxInfo[c->tid].curFullCount, auxInfo[c->tid].fp_fullMatch);
        }
        if (auxInfo[c->tid].curOtherCount > 0)
        {
                if (auxInfo[c->tid].fp_others == NULL)
                {
                        auxInfo[c->tid].fp_others = fopen(auxInfo[c->tid].otherFileName, "w");
                        assert(auxInfo[c->tid].fp_others != NULL);
                }
                fwrite(auxInfo[c->tid].Or, sizeof(otherRec), auxInfo[c->tid].curOtherCount, auxInfo[c->tid].fp_others);
        }

	bam_hdr_destroy(h);
	bam_destroy1(b);
}
int main(int argc, char **argv)
{
	int i=0;
	if (argc != 4)
	{
		printf ("Usage: a.out <numThreads> <bamFile> <OutPrefix>\n");
		return 1;
	}
	int numThreads = atoi(argv[1]);
	char *bamFileName = argv[2];
	char *outPrefix = argv[3];
	char fName[500];
	htsFile *fp1=hts_open(bamFileName, "r");
	bam_hdr_t *hdr1 = sam_hdr_read(fp1);
	auxInfo = (contigInfo *)malloc(hdr1->n_targets * sizeof(contigInfo));
        for (i=0; i<hdr1->n_targets; i++)
        {
                sprintf (fName, "%s_%s_0.aeb", outPrefix, hdr1->target_name[i]);
                auxInfo[i].fullFileName = (char *)malloc(strlen(fName)+1);
                strcpy(auxInfo[i].fullFileName, fName);
                sprintf (fName, "%s_%s_0.aib", outPrefix, hdr1->target_name[i]);
                auxInfo[i].otherFileName = (char *)malloc(strlen(fName)+1);
                strcpy(auxInfo[i].otherFileName, fName);
                auxInfo[i].fp_fullMatch = NULL;
                auxInfo[i].fp_others = NULL;
                auxInfo[i].curFullCount=0;
                auxInfo[i].curOtherCount=0;
                auxInfo[i].outputChunkSize=100000;
                auxInfo[i].Fr = (fullRec *)malloc(auxInfo[i].outputChunkSize*sizeof(fullRec));
                auxInfo[i].Or = (otherRec *)malloc(auxInfo[i].outputChunkSize*sizeof(otherRec));
                assert (auxInfo[i].Fr != NULL);
                assert (auxInfo[i].Or != NULL);
        }
	
	#pragma omp parallel for private(i) num_threads(numThreads)
	for (i=0; i<hdr1->n_targets; i++)
	{
		htsFile *fp=hts_open(bamFileName, "r");
		assert (fp!=NULL);
		hts_idx_t *idx = sam_index_load(fp, bamFileName);  // load the index
	        if (idx == NULL) {
	            printf("can't load index for \"%s\"", argv[1]);
	        }
		bam_hdr_t *hdr = sam_hdr_read(fp);
		bam1_t *b = (bam1_t*)calloc(1, sizeof(bam1_t));

		hts_itr_t *iter = sam_itr_queryi(idx, i, 1, hdr->target_len[i]);
		if (iter == NULL) {
                	printf("can't parse region \"%s\"\n", hdr->target_name[i]);
            	}
		int ret = sam_itr_next(fp, iter, b);
		if (ret == -1) continue;
		printf ("Contig - %s\n",hdr->target_name[i]);
		readBamRec(b, hdr, iter, fp, auxInfo[i].fp_fullMatch, auxInfo[i].fp_others);
		if (auxInfo[i].fp_fullMatch != NULL)
			fclose (auxInfo[i].fp_fullMatch);
		if (auxInfo[i].fp_others != NULL)
			fclose (auxInfo[i].fp_others);
		iter=NULL;
		hts_idx_destroy(idx);
		hts_close (fp);
	}
	hts_close (fp1);
	
	return 0;
}
