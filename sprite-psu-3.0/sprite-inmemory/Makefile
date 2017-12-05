include Make.inc
all: sprite-inmemory

sprite-inmemory: parsnip.o sprite.o
	$(CC) $(CFLAGS) -o prune-genreadidx ../prune/prune-genreadidx.c -lz
ifeq ($(MPI),1)
	$(MPICC) $(CFLAGS) -DUSE_MPI parsnip.o $(BWA_DIR)/kopen.o sprite.o -o sprite-inmemory -L$(BWA_DIR) -lbwa -lm -lz -lpthread
else
	$(CC) $(CFLAGS) parsnip.o $(BWA_DIR)/kopen.o sprite.o -o sprite-inmemory -L$(BWA_DIR) -lbwa -lm -lz -lpthread
endif
sprite.o: sprite.c
ifeq ($(MPI),1)
	$(MPICC) $(CFLAGS) -DUSE_MPI -c -I$(BWA_DIR) sprite.c
else
	$(CC) $(CFLAGS) -c -I$(BWA_DIR) sprite.c
endif
parsnip.o: parsnip.c
ifeq ($(MPI),1)
	$(MPICC) $(CFLAGS) -DUSE_MPI -c parsnip.c
else
	$(CC) $(CFLAGS) -c parsnip.c
endif
clean:
	rm -f prune-genreadidx parsnip.o sprite.o sprite-inmemory
