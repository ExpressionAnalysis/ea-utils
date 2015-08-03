# Sam-stats statistics by name #

### reads ###

Incremented for every sequence-containing line in the sam file, regardless of whether it represents an alignment.   for some files, this is not actually the number of reads.  indeed, this may be a poor name for this stat

### version ###

The version of sam-stats used.

### mapped reads ###

If duplicate tracking was enabled via -D, then this attempts to recapitulate the number of unique, mapped, probe-id's in the original sam file.  It is multiplied by 2 for paired-end data with duplicate read id's.  The idea is that if you divide this by the number of reads in the fastq you aligned (possibly from the output of fastq-stats), you will get an accurate "percentage of reads aligned" statistic.

If duplicate tracking is not enabled, then the result is the number of mapped reads in the sam file.

The definition of "mapped" is something with a non-negative position, and a "non-asterisk" cigar string.

### pct ambiguous ###

Formula is:

100 **<number of probe id's aligning more than once> / <number of probe id's>**

"Aligning more than once" means more than 2 alignments for paired end and more than 1 alignment for single-end.

_This does NOT, currently, use the BWA XA tag, and support for this will be added in the near future._

### max dup align ###

The maximum number of alignments for a single probe-id.  IE: if seq4455 aligns 50 times, and this is more than any other sequence, then the result would be 50.   For paired-end reads one is subtracted, since you should get two alignments with the same probe-id for each read, but anything more than that is considered "ambiguous".

One could argue, rightly, that it should be divided by two for paired end.  If this is changed, it would require a major version upgrade, or rename of the statistic since it would invalidate historical comparisons.


### singleton mappings ###

Number of reads that aligned only once, without a pair.  Only tracked if -D is enabled.

### total mappings ###

The number of lines in the sam file that were "mapped" reads.

### mapped bases ###

The total number of bases, derived from the length of the sequence strings, in the sam file that were "mapped" reads.

### library:paired-end ###

And indicator that "paired end" semantics were used when running sam-stats, because there were alignments that exhibited these characteristics.

### discordant mates ###
The number of read-pairs that were aligned to different chromosomes.

### distant mates ###
The number of read-pairs that were aligned more than 50k bases apart

### phred ###
The phred scale used when analyzing this sam file.  Should always be 33, unless you're looking at some old Illumina aligned output.


### forward ###
The number of lines in the sam file that were aligned to the "forward" strand.  No accounting is done on duplicates.

### reverse ###
The number of lines in the sam file that were aligned to the "reverse" strand.  No accounting is done on duplicates.

### len (max/min/mean) ###
"len max" is always output.   "len mean, len stdev" are only output if the length of the sequence reads vary.  Calculations are based on the the length of the (possibly hard-clipped) sequence in the sam file.

### mapq (mean/stdev/Q1/median/Q3) ###

These are all statistics on mapped reads only, with no attempt to de-duplicate multiply aligned reads before calculating the statistic.   As a result, a single read that is represented 1000 times in the sam results can skew the mapping quality significantly lower (since they will all be zeroes).

The statistics presented are: mean and standard deviation, (Q1) first-quartile (R-compatible algorithm used), median, and (Q3) upper-quartile.

### snp rate ###
The total number of mismatch bases, divided by the total number of mapped based.  Calculated using NM tags.

### ins rate ###
The total number of inserted bases, divided by the total number of mapped based.  Calculated using cigar strings.

### del rate ###
The total number of deleted bases, divided by the total number of mapped based.  Calculated using cigar strings.

### insert mean/stdev/Q1/median/Q3 ###
These are the "isize" statistics.   The mean and standard deviation are ["trimmed" by 10%](http://en.wikipedia.org/wiki/Truncated_mean), in accordance with best-practices for producing meaningful fragment size distributions.

You will have to add the mean length of read2 to this value in order to get the total fragment size for programs like Mosaik and RSEM that require it, or subtract the mean read1 size for programs like tophat that require the "between fragment" size (even if negative).

Quartiles are computed with an R-compatible algorithm

### base qual mean/stdev ###

Total quality score of every mapped base, divided by the total number of mapped bases.

### %A/%T/%C/%G/%N ###

Total number of mapped bases of each type, divided by the total number of mapped bases.

### num ref seqs ###

The number of reference sequences in the original SAM headers.

### num ref aligned ###

The number of unique reference sequences that have at least one alignment.

### median skew ###

Only enabled in "rna mode", the median skewness of all transcripts with coverage.  An "overall bias" metric, useful for detecting 3'bias.

### median coverage ###

Median of nonzero coverage (see rna mode below)

### median coverage cv ###

Median of nonzero coverage variability (see rna mode below)

### %chr1, %gene, and signature values ###

For each mapped reference chromosome or gene, the total number of mapped bases is totaled.   The percentages are the total mapped bases for that gene or chromosome divided by the total mapped bases.

If "-A" is used, up to a million reference sequences are tracked and output.

In addition to outputting a percentage, and only if the sam file had a header, a signature representing a histogram of alignments is output.

Each chromosome/gene is divided in to 30 evenly sized "buckets".

The values in the signature are the log2 of the # of mapped reads in  each position bucket times the length of those reads + ascii '0'.

Long reads caveat: There is no attempt to break-up long reads that align across multiple regions however.  So these signatures are not (yet) useful for very long reads.   If a read's first position maps to histogram bucket 1, but spans across to bucket 2, all of the bases are currently counted in bucket 1.   For short reads, this is not significant.


### "RNA-seq" mode coverage matrix ###

A "coverage" file is output with one row per reference sequence, the _7 columns_ are: transcript id, length, count, coverage percentage, skewness (bias), and coverage cv, and a signature (see above).

  * Coverage: statistics are 1x, and approximate, and are intended to measure exon bias - not appropriate for assembly!

  * Skewness: is the standard pearsons 3rd order skewness applies to the alignment positions across the transcript.

  * Coverage CV: the coefficient of variation of the coverage histogram (can be computed from the signature).