# Usage #

```
Usage: fastq-join [options] <read1.fq> <read2.fq> [mate.fq] -o <read.%.fq>

Joins two paired-end reads on the overlapping ends.

Options:

-o FIL          See 'Output' below
-v C            Verifies that the 2 files probe id's match up to char C
                  use '/' for Illumina reads
-p N            N-percent maximum difference (8)
-m N            N-minimum overlap (6)
-r FIL          Verbose stitch length report
-R              No reverse complement
-V              Show version

Output:

  You can supply 3 -o arguments, for un1, un2, join files, or one
argument as a file name template.  The suffix 'un1, un2, or join' is
appended to the file, or they replace a %-character if present.

  If a 'mate' input file is present (barcode read), then the files
'un3' and 'join2' are also created.

  Files named ".gz" are assumed to be compressed, and can be 
read/written as long as "gzip" is in the path.
```

## Etc ##

This uses our sqr(distance)/len for anchored alignment quality algorithm.  It's a good measure of anchored alignment quality, akin (in my mind) to squared-deviation for means.

## Overlapping Bases ##

### When the bases match ###

The higher quality base is used, and it is increased by up to 3

### When the bases don't match ###

If one quality is greater than "3" (50%), then the the resulting quality is the difference between the two qualities (reduced quality due to mismatch), or "3" )(50%), whichever is greater.

#### Examples: ####

```
40 vs 3 = 37 : second base has low quality... doesn't change top by much
40 vs 40 = 3 : two equal quality bases that don't match = qual of 3
2 vs 2 = 2 : neither base has a high quality
```

#### Some caveats: ####

Illumina's quality scores are not accurate and estimates vary by chemistry and sequencer.   I would recommend using a profiling tool,  on PhiX, and adjusting your qualities using the results of the tool.

For example the quality score "2" has a true quality, typically of "11" ... this is Illumina's code for "quality estimation failure".   The quality scores at the high end ("34-40") are often overestimates.