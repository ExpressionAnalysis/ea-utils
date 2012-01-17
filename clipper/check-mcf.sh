#!/bin/bash

g++ fastq-mcf.c fastq-lib.cpp -o fastq-mcf.ex

ok=/opt/bin/fastq-mcf
new=./fastq-mcf.ex
err=0

err()
{
	err=1	
} 

if [ "$1" = "-b" ]; then
	list="new ok"
else
	list="new"
fi

trap err ERR

set -o xtrace
for v in $list; do
	eval prog=\$$v
	${prog} test.fa test1.fq > test1.$v.out 2> test1.$v.err || echo err during 1 $v
	${prog} test.fa test2.fq > test2.$v.out 2> test2.$v.err || echo err during 2 $v
	${prog} test.fa test1.fq -o test3.$v.out > test3.$v.err || echo err during 3 $v
	${prog} -L 72 -f test.fa test4.fq1 test4.fq2 -o test4.$v.out -o test4.$v.out2 > /dev/null || echo err during 4 $v
# check skip saving
	${prog} test.fa test1.fq -S -o test5.$v.out > test5.$v.err || echo err during 5 $v
# check gzipping
	${prog} test.fa test1.fq -o test6.$v.out.gz > test6.$v.err || echo err during 6 $v
	mv test5.$v.out.skip test5.$v.out
done
set +o xtrace	

for n in test1 test2 test3 test4 test5; do
	echo $n
	diff $n.new.out $n.ok.out
	[[ -e $n.ok.err ]] && diff $n.new.err $n.ok.err
done

n=test6
echo $n
zdiff $n.new.out.gz $n.ok.out.gz

shopt -s extglob

rm test?.@(new).@(out|err)

if [ $err == 1 ]; then
	echo NOT OK, some differences
	exit 1
else
	echo OK, all tests passed
	exit 0
fi
