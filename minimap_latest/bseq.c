//#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "bseq.h"
#include "kseq.h"
KSEQ_INIT(FILE *, fread)

extern unsigned char seq_nt4_table[256];

struct bseq_file_s {
	int is_eof;
	FILE *fp;
	kseq_t *ks;
};

bseq_file_t *bseq_open(const char *fn)
{
	bseq_file_t *fp;
	FILE * f;
	f = fn && strcmp(fn, "-")? fopen(fn, "r") : fdopen(fileno(stdin), "r");
	if (f == 0) return 0;
	fp = (bseq_file_t*)calloc(1, sizeof(bseq_file_t));
	fp->fp = f;
	fp->ks = kseq_init(fp->fp);
//	printf ("ofs %lu\n", ftell(f));
	return fp;
}

void bseq_close(bseq_file_t *fp)
{
	kseq_destroy(fp->ks);
	fclose(fp->fp);
	free(fp);
}

void bseq_seek(bseq_file_t *fp, size_t ofs)
{
	fseek (fp->fp, ofs, SEEK_SET);
//	printf ("ofs %lu\n", ftell(fp->fp));
}

int bseq_tell(bseq_file_t *fp)
{
	return (int)(ftell(fp->fp));
}
bseq1_t *bulk_read(bseq_file_t *fp, int chunk_size, int *n_, int *sizeToRead, char **buf, size_t *bufSize)
{
	if (*bufSize < chunk_size)
	{
		*buf = (char *)realloc (*buf, chunk_size);
		assert (buf != NULL);
		*bufSize = chunk_size;
	}
//	printf ("chunk_size %d\n", chunk_size);
	int size = 0, m, n;
	bseq1_t *seqs;
	kseq_t *ks = fp->ks;
	m = n = 0; seqs = 0;
//	char *tmpBf = (char *)malloc (chunk_size);
//	assert (tmpBf != NULL);
//	printf ("bulk_read tmpBf %p size %d\n", tmpBf, chunk_size);
	char *tmpBf = *buf;
	assert (fread (tmpBf, 1, chunk_size, fp->fp) == chunk_size);
	int i;
	for (i=0; i<chunk_size;)
	{
		assert (tmpBf[i] == '@');
		i++;

		bseq1_t *s;
		if (n >= m) {
			m = m? m<<1 : 256;
//			printf ("bulk_read seqs %p m %lu bytes\n", seqs, m * sizeof(bseq1_t));
			seqs = (bseq1_t*)realloc(seqs, m * sizeof(bseq1_t));
		}
		s = &seqs[n];

		
		s->name = &tmpBf[i];
		s->l_name=0;
		while (!isspace (tmpBf[i])) 
		{
			s->l_name++;
			i++;
		}
		if (tmpBf[i] != '\n')
		{
			tmpBf[i] = '\0';
			i++;
			while (tmpBf[i] != '\n') i++;
		}
		assert (tmpBf[i] == '\n');
		tmpBf[i] = '\0';
		i++;
		//End Of name line
		s->l_seq=0;
		s->seq = &tmpBf[i];
		while (tmpBf[i] != '\n')
		{
			s->l_seq++;
			i++;
		}
		tmpBf[i] = '\0';
		i++;
		//End of SEQ line
		while (tmpBf[i] != '\n')i++;
		i++;
		//End of + line
		s->qual = &tmpBf[i];
		i+=s->l_seq;
		assert (tmpBf[i] == '\n');//End of Read
		tmpBf[i] = '\0';
		i++;//Point to next read
//		printf ("Read: name(%d) *%s*\n seq(%d) *%s*\n qual *%s*\n", s->l_name, s->name, s->l_seq, s->seq, s->qual);
		
		n++;
	}
	
	if (sizeToRead != NULL)
		*sizeToRead = 0;
	if (n == 0) fp->is_eof = 1;
	*n_ = n;
	return seqs;
}
bseq1_t *bseq_read(bseq_file_t *fp, int chunk_size, int *n_, int *sizeToRead)
{
//	printf ("chunk_size %d\n", chunk_size);
	int size = 0, m, n;
	bseq1_t *seqs;
	kseq_t *ks = fp->ks;
	m = n = 0; seqs = 0;
	char *tmpBf = (char *)malloc (chunk_size);
	printf ("bseq_read tmpBf %p size %d\n", tmpBf, chunk_size);
	assert (tmpBf != NULL);
	fread (tmpBf, 1, chunk_size, fp->fp);
	
/*	while ((chunk_size > 0) && (kseq_read(ks) >= 0)) {
		bseq1_t *s;
		assert(ks->seq.l <= INT32_MAX);
		if (n >= m) {
			m = m? m<<1 : 256;
			seqs = (bseq1_t*)realloc(seqs, m * sizeof(bseq1_t));
		}
		s = &seqs[n];
//		s->name = strdup(ks->name.s);
//		s->seq = strdup(ks->seq.s);
		s->l_seq = ks->seq.l;
		s->l_name = ks->name.l;
		if (sizeToRead != NULL)
		{
//	                s->qual = strdup(ks->qual.s);//vasu
			size += ks->name.l + ks->comment.l + ks->plus_len + 2*seqs[n].l_seq + 4 + 2 + (ks->comment.l > 0);
		}
		else
			size += seqs[n].l_seq;
		n++;
//		printf ("%d - *%s* *%s* - %d, %d, %d, %d\n", size, ks->name.s, ks->comment.s, ks->name.l, ks->comment.l, ks->plus_len, ks->seq.l);
		if (size >= chunk_size) break;
	}*/
	seqs = (bseq1_t *)tmpBf;
	if (sizeToRead != NULL)
		*sizeToRead -= size;
	if (n == 0) fp->is_eof = 1;
	*n_ = n;
	return seqs;
}
bseq1_t *bseq_read2(bseq_file_t *fp, int n_)
{
	int size = 0, m, n;
	bseq1_t *seqs;
	kseq_t *ks = fp->ks;
	m = n = 0; seqs = 0;
	while (kseq_read(ks) >= 0) {
		bseq1_t *s;
		assert(ks->seq.l <= INT32_MAX);
		if (n >= m) {
			m = m? m<<1 : 256;
			seqs = (bseq1_t*)realloc(seqs, m * sizeof(bseq1_t));
		}
		s = &seqs[n];
		s->name = strdup(ks->name.s);
		s->seq = strdup(ks->seq.s);
                s->qual = strdup(ks->qual.s);//vasu
		s->l_seq = ks->seq.l;
		s->l_name = ks->name.l;
		size += ks->name.l + ks->comment.l + ks->plus_len + 2*seqs[n].l_seq + 4 + 2 + (ks->comment.l > 0);
		n++;
//		printf ("%d - *%s* *%s* - %d, %d, %d, %d\n", size, ks->name.s, ks->comment.s, ks->name.l, ks->comment.l, ks->plus_len, ks->seq.l);
		if (n >= n_) break;
	}
	if (n == 0) fp->is_eof = 1;
	assert (n_ == n);
	return seqs;
}

int bseq_eof(bseq_file_t *fp)
{
	return fp->is_eof;
}
