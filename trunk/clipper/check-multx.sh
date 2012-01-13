#!/bin/bash

g++ fastq-multx.c && ./a.out -v ' ' -l /mnt/scratch/all_barcodes.txt mxtest_2.fastq mxtest_1.fastq mxtest_3.fastq -o n/a -o mxout_%_1.fq.tmp -o mxout_%_2.fq.tmp > t1

g++ fastq-multx.c && ./a.out -v ' ' -l /mnt/scratch/all_barcodes.txt mxtest_2.fastq mxtest_1.fastq mxtest_3.fastq -o n/a -o mxout_%_1.fq.tmp.gz -o mxout_%_2.fq.tmp.gz > t2

for x in LB2 LB4 LB5 LB6; do
	if ( ! diff  mxout_${x}_1.fq.tmp  mxout_${x}_1.fq > /dev/null ); then
		echo FAIL: diff  mxout_${x}_1.fq.tmp  mxout_${x}_1.fq
		exit 1
	fi
done

for x in t1 t2; do
if ( ! diff  mxout.$x.txt $x ); then
	echo FAIL: diff t1 mxout.txt
	exit 1
fi
done

rm mx*.tmp mx*.tmp.gz t1 t2

echo OK
