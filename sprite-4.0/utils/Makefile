include ../Make.inc
all: aebaib2txt bam2aebaib

aebaib2txt: aebaib2txt.c
	$(CC) $(CFLAGS) -o aebaib2txt aebaib2txt.c 
bam2aebaib: bam2aebaib.c
	${CC} -I$(HTSLIB_DIR)/htslib $(CFLAGS) $(OPENMP_FLAG)  -o bam2aebaib bam2aebaib.c -L$(HTSLIB_DIR)  -lhts -lz
clean:
	rm -f aebaib2txt bam2aebaib
