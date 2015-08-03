# Introduction #

Tool for computing statistics from (possibly compressed) SAM or BAM files.

See SamStatsDetails for more info on the output fields.

# Usage #

```
Usage: sam-stats [options] [file1] [file2...filen]
Version: 1.32

Produces lots of easily digested statistics for the files listed

Options (default in parens):

-D             Keep track of multiple alignments (slower!)
-M             Only overwrite if newer (requires -x, or multiple files)
-A             Report all chr sigs, even if there are more than 1000
-R FIL         RNA-Seq stats output (coverage, 3' bias, etc)
-B             Input is bam, dont bother looking at magic
-x FIL         File extension for multiple files (stats)
-b INT         Number of reads to sample for per-base stats (1M)
-S INT         Size of ascii-signature (30)
-z             Don't fail when zero entries in sam

OUTPUT:

If one file is specified, then the output is to standard out.  If
multiple files are specified, or if the -x option is supplied,
the output file is <filename>.<ext>.  Default extension is 'stats'.

Complete Stats:

  <STATS>           : mean, max, stdev, median, Q1 (25 percentile), Q3
  reads             : # of entries in the sam file, might not be # reads
  phred             : phred scale used
  bsize             : # reads used for qual stats
  mapped reads      : number of aligned reads (unique probe id sequences)
  mapped bases      : total of the lengths of the aligned reads
  forward           : number of forward-aligned reads
  reverse           : number of reverse-aligned reads
  secondary         : # of entries with 256-bit set
  snp rate          : mismatched bases / total bases
  ins rate          : insert bases / total bases
  del rate          : deleted bases / total bases
  pct mismatch      : percent of reads that have mismatches
  len <STATS>       : read length stats, ignored if fixed-length
  mapq <STATS>      : stats for mapping qualities
  insert <STATS>    : stats for insert sizes
  %<CHR>            : percentage of mapped bases per chr, followed by a signature

Subsampled stats (1M reads max):
  base qual <STATS> : stats for base qualities
  %A,%T,%C,%G       : base percentages

Meaning of the per-chromosome signature:
  A ascii-histogram of mapped reads by chromosome position.
  It is only output if the original SAM/BAM has a header. The values
  are the log2 of the # of mapped reads at each position + ascii '0'.
```

DETAILED OVERVIEW OF STATISTICS:

Click SamStatsDetails for more information on each stat, how it's calculated and what it means.