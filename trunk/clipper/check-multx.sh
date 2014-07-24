#!/bin/bash

# TODO: THIS IS AWFUL... MAKE IT REAL

set -o xtrace
g++ fastq-multx.cpp fastq-lib.cpp


./a.out -v ' ' -l master-barcodes.txt mxtest_2.fastq mxtest_1.fastq mxtest_3.fastq -o n/a -o mxout_%_1.fq.tmp -o mxout_%_2.fq.tmp > t1 2> e1
./a.out -v ' ' -l master-barcodes.txt mxtest_2.fastq mxtest_1.fastq mxtest_3.fastq -o n/a -o mxout_%_1.fq.tmp.gz -o mxout_%_2.fq.tmp.gz > t2 2> e2
./a.out -g mxtest_2.fastq mxtest_1.fastq mxtest_3.fastq -o mzout_%_1.fq.tmp -o mzout_%_2.fq.tmp > t3 2> e3
./a.out -v ' ' -l master-barcodes.txt mxtest-h_1.fastq mxtest-h_2.fastq -o n/a -o mwout_%_1.fq.tmp -o mwout_%_2.fq.tmp > t4 2> e4

for x in LB2 LB4 LB5 LB6; do
	if ( ! diff  mxout_${x}_1.fq.tmp  mxout_${x}_1.fq > /dev/null ); then
		echo FAIL: diff  mxout_${x}_1.fq.tmp  mxout_${x}_1.fq
		err=1
	fi
done

for x in t1 t2 t3 t4; do
if ( ! diff  mxout.$x.txt $x ); then
	echo FAIL: diff $x mxout.$x.txt
    err=1
fi
done

if [ $err ]; then
    exit 1
fi

rm mw*.tmp mz*.tmp mx*.tmp mx*.tmp.gz t1 t2 t3 e1 e2 e3

echo OK
