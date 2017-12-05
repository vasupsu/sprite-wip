include ../Make.inc
BWA_DIR=../bwa-0.7.12

all: prune-genreadidx prune-bwamem bam2aebaib

prune-genreadidx: prune-genreadidx.c
	$(CC) $(CFLAGS) -o prune-genreadidx prune-genreadidx.c -lz

prune-bwamem: faiTrie.o prune-bwamem.o
ifeq ($(MPI),1)
	$(MPICC) $(CFLAGS) faiTrie.o prune-bwamem.o $(BWA_DIR)/kopen.o -o prune-bwamem -L$(BWA_DIR) -lbwa -lm -lz -lpthread
else
	$(CC) $(CFLAGS) faiTrie.o prune-bwamem.o $(BWA_DIR)/kopen.o -o prune-bwamem -L$(BWA_DIR) -lbwa -lm -lz -lpthread
endif

faiTrie.o: faiTrie.c
	$(CC) $(CFLAGS) -c faiTrie.c

prune-bwamem.o: prune-bwamem.c
ifeq ($(MPI),1)
	$(MPICC) $(CFLAGS) -DUSE_MPI -c -I$(BWA_DIR) prune-bwamem.c
else
	$(CC) $(CFLAGS) -c -I$(BWA_DIR) prune-bwamem.c
endif

bam2aebaib: bam2aebaib.c
	$(CC) -I$(HTSLIB_DIR)/htslib $(CFLAGS) $(OPENMP_FLAG)  -o bam2aebaib bam2aebaib.c -L$(HTSLIB_DIR)  -lhts -lz

clean:
	rm -f *.o prune-genreadidx prune-bwamem bam2aebaib
