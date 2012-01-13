#!/bin/bash -e

g++ fastq-clipper.c

ok=/opt/bin/fastq-clipper
new=./a.out

# comparre
fastx_clipper -i test1.fq -a AGTCCCGTAC -o test1.fx.out

for v in new ok; do
	eval prog=\$$v
	${prog} test1.fq AGTCCCGTAC > test1.$v.out 2> test1.$v.err
	${prog} test2.fq AGTCCCGTAC > test2.$v.out 2> test2.$v.err
	${prog} test1.fq AGTCCCGTAC -o test3.$v.out > test3.$v.err
	diff test1.$v.out test1.fx.out > test4.$v.out && true
done
	

for n in test1 test2 test3 test4; do
	echo $n
	diff $n.new.out $n.ok.out
	[[ -e $n.ok.err ]] && diff $n.new.err $n.ok.err
done

shopt -s extglob

rm test?.@(new|ok).@(out|err)

echo OK, all tests passed
