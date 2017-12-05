CC=gcc
MPICC=mpicc
# MPI=1 to compile for multi-node run using MPI
MPI=1
CFLAGS=-g -Wall -O2
OPENMP_FLAG=-fopenmp
# specify path to BWA (for fastmap.o and libbwa.a)
BWA_DIR=bwa-0.7.12
# specify path to htslib
HTSLIB_DIR=/work/03201/tg824998/samtools-1.1/htslib-1.1
