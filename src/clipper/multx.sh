#!/bin/sh
#g++ -g fastq-mcf.c -o fastq-mcf.ex && ./fastq-mcf.ex $*
g++ -g fastq-multx.c -o fastq-multx.ex && gdb --eval-command=run --args ./fastq-multx.ex $*
