#!/bin/sh

g++ fastq-mcf.c -o fastq-mcf.ex

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
	${prog} test.fa test1.fq > test1.$v.out 2> test1.$v.err
	${prog} test.fa test2.fq > test2.$v.out 2> test2.$v.err
	${prog} test.fa test1.fq -o test3.$v.out > test3.$v.err
	${prog} -L 72 -f test.fa test4.fq1 test4.fq2 -o test4.$v.out -o test4.$v.out2
done
set +o xtrace	

for n in test1 test2 test3 test4; do
	echo $n
	diff $n.new.out $n.ok.out
	[[ -e $n.ok.err ]] && diff $n.new.err $n.ok.err
done

shopt -s extglob

rm test?.@(new).@(out|err)

if [ $err == 1 ]; then
	echo NOT OK, some differences
	exit 1
else
	echo OK, all tests passed
	exit 0
fi
