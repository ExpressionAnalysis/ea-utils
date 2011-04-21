#!/bin/sh -e

path=`cat install.path`

g++ -O3 fastq-clipper.c -o $path/fastq-clipper&
g++ -O3 fastq-mcf.c -o $path/fastq-mcf&
g++ -O3 fastq-multx.c -o $path/fastq-multx&
g++ -O3 fastq-join.c -o $path/fastq-join&
wait

echo OK
