#!/bin/bash
BWA_exec=../bwa-0.7.12/bwa
SAMTOOLS_exec=samtools/samtools

echo " ****** Getting some test data ..."
wget https://psu.box.com/shared/static/edi1kjukyeku2zhf8yaelec5pb4t40id.fa
mv edi1kjukyeku2zhf8yaelec5pb4t40id.fa ref.fa
wget https://psu.box.com/shared/static/3505dplglwcg4dx3xecwjk1zlrmdwv53.fq
mv 3505dplglwcg4dx3xecwjk1zlrmdwv53.fq reads_1.fa
wget https://psu.box.com/shared/static/ugln06us9kubrdkdps1dkbn9njkx0ov8.fq
mv ugln06us9kubrdkdps1dkbn9njkx0ov8.fq reads_2.fa
echo " ****** Done downloading test data."

echo " ****** Building index using BWA ..."
$BWA_exec index ref.fa
$SAMTOOLS_exec faidx ref.fa
echo " ****** Done building index."

echo " ****** In-Memory SPRITE(1 thread) ..."
mkdir outdir
../prune-genreadidx reads_1.fa 100
../sprite-inmemory -t 1 -o outdir/out ref.fa reads_1.fa reads_2.fa
echo " ****** Done running In-memory SPRITE."

#echo " ****** Removing intermediate files"
rm -f outdir/out_* ref.fa.* reads_1.fa.idx
#echo " ****** Done cleaning up."
