include ../Make.inc
all: sampa parsnip

sampa: sampa.c
ifeq ($(MPI),1)
	$(MPICC) $(CFLAGS) -DUSE_MPI -o sampa sampa.c
else
	$(CC) $(CFLAGS) -o sampa sampa.c
endif
parsnip: parsnip.c
ifeq ($(MPI),1) 
	$(MPICC) $(CFLAGS) -DUSE_MPI -o parsnip parsnip.c -lm -lpthread
else
	$(CC) $(CFLAGS) -o parsnip parsnip.c -lm -lpthread
endif
clean:
	rm -f sampa parsnip
