# Introduction #

The idea behind this is to reduce the amount of "piping" going on in a pipeline.   A lot of time, disk space and nail-chewing is spent keeping files in sync, figuring out what barcodes are on what samples, etc.  The goal of this program is to make it easier to demultiplex possibly paired-end sequences, and also to allow the "guessing" of barcode sets based on master lists of barcoding protocols (fluidigm, truseq, etc.)

# Usage #

```
Usage: fastq-multx [-g|-l] <barcodes.fil> <read1.fq> -o r1.%.fq [mate.fq -o r2.%.fq] ...

Output files must contain a '%' sign which is replaced with the barcode id in the barcodes file.

Barcodes file looks like this:

<id1> <sequence1>
<id2> <sequence2> ...

Default is to guess the -bol or -eol based on clear stats.

If -g is used, then it's parameter is an index lane, and frequently occuring sequences are used.

If -l is used then all barcodes in the file are tried, and the *group* with the *most* matches is chosen.

Grouped barcodes file looks like this:

<id1> <sequence1> <group1>
<id1> <sequence1> <group1>
<id2> <sequence2> <group2>...

Mated reads, if supplied, are kept in-sync

Options:

-o FIL1 [FIL2]  Output files (one per input, required)
-g FIL          Determine barcodes from indexed read FIL
-l FIL          Determine barcodes from any read, using FIL as a master list
-b              Force beginning of line
-e              Force end of line
-x              Don't trim barcodes before writing
-n              Don't execute, just print likely barcode list
-v C            Verify that mated id's match up to character C ('/' for illumina)
-m N            Allow up to N mismatches, as long as they are unique
```


Files named ".gz" are assumed to be compressed, and can be
read/written as long as "gzip" is in the path.

# Example 1 #

# this example will read/output files that are gzipped, since -B is used and only 1 sequence files is present, it will look for barcodes on the "ends" of the sequence and will tell you which end it found them on

`fastq-multx -B barcodes.fil seq.fastq.gz -o %.fq.gz`

Contents of barcodes.fil:

```
mock_a ACCC
salt_a CGTA
mock_b  GAGT
salty_b  CGGT
```

# Example 2 #

# this example will first determine which "barcode group" to use, will select the most likely set of barcodes from that file, and will then proceed as if only that set was specified.   this allows for a single pipeline that works with multiple technologies

`fastq-multx -l barcodes.grp seq2.fastq.gz seq1.fastq.gz -o n/a -o out%.fq`

Contents of barcodes.grp:

```
id      seq     style
LB1     ATCACG  TruSeq
LB2     CGATGT  TruSeq
LB3     TTAGGC  TruSeq
LB4     TGACCA  TruSeq
LB5     ACAGTG  TruSeq
A01_01  TAGCTTGT        Fluidigm
B01_02  CGATGTTT        Fluidigm
C01_03  GCCAATGT        Fluidigm
D01_04  ACAGTGGT        Fluidigm
E01_05  ATCACGTT        Fluidigm
```

Standard error will output:

`Using Barcode Group: TruSeq on File: seq2.fastq.gz (start), Threshold 0.59%`

This indicated that The LB1-LB5 barcodes will be used, and that the filess will be named LB1-LB5, and that the barcode was at the "start" of the reads in the seq2 file.


Example of Nextera/Dual-Indexed input:

```
id      seq     style
D708_508        TAATGCGC-GTACTGAC       TruSeq RNA
D709_501        CGGCTATG-TATAGCCT       TruSeq RNA
D709_502        CGGCTATG-ATAGAGGC       TruSeq RNA
D709_503        CGGCTATG-CCTATCCT       TruSeq RNA
D709_504        CGGCTATG-GGCTCTGA       TruSeq RNA
D709_505        CGGCTATG-AGGCGAAG       TruSeq RNA
D709_506        CGGCTATG-TAATCTTA       TruSeq RNA
D709_507        CGGCTATG-CAGGACGT       TruSeq RNA
D709_508        CGGCTATG-GTACTGAC       TruSeq RNA
D710_501        TCCGCGAA-TATAGCCT       TruSeq RNA
D710_502        TCCGCGAA-ATAGAGGC       TruSeq RNA
D710_503        TCCGCGAA-CCTATCCT       TruSeq RNA
D710_504        TCCGCGAA-GGCTCTGA       TruSeq RNA
D710_505        TCCGCGAA-AGGCGAAG       TruSeq RNA
D710_506        TCCGCGAA-TAATCTTA       TruSeq RNA
D710_507        TCCGCGAA-CAGGACGT       TruSeq RNA
```