#!/bin/bash

./fastq-join -p 20 -m 5 test-overlap/test_1x test-overlap/test_2x -o test-overlap/a.testx.

./fastq-join -p 20 -m 5 test-overlap/test_1 test-overlap/test_2 -o test-overlap/a.test.

./fastq-join -p 20 -m 5 test-overlap/test_1x test-overlap/test_2x -o test-overlap/b.testx. -x

./fastq-join -p 20 -m 5 test-overlap/test_1 test-overlap/test_2 -o test-overlap/b.test. -x

