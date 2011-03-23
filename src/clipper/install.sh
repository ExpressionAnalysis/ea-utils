#!/bin/sh -e

g++ -O3 fastq-clipper.c -o /opt/bin/fastq-clipper
g++ -O3 fastq-mcf.c -o /opt/bin/fastq-mcf
g++ -O3 fastq-multx.c -o /opt/bin/fastq-multx
g++ -O3 fastq-join.c -o /opt/bin/fastq-join
