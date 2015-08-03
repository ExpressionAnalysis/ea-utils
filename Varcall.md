# Introduction #

varcall does two tasks:

  * Calculates error rate statistics in a pileup
  * Outputs variants per position


The tool is not complete.

In particular: the VCF output is somewhat broken.   The native output formats now work fine, and the stats/outputs collection and use are also working well.

One of the advantages of using varcall over something like GATK is the extensive parameterization, allowing varcall to be used in any situation, this is also a shortcoming... since users must know the meanings of every parameter.  In addition, varcall seems to be particularly fast.

Recently we added variation specific error rate collection, and it is used, instead of a the global error rate, in the calculation of pvals.   This vastly reduces false positives, etc.

By the time we reach version 1, all of these should be fixed, including better docs, and better "auto-configuration" based on sampling.  If you are using varcall before then, carefully experiment with parameters and compare to other tools.

Version 1 should also have the ability to do a stats collection across a group of samples, and then subsequent calling using those stats.   The advantage here is that for some amplicon/pcr assays, there are insufficient locii to get an accurate error variation estimate per-sample.


# Usage #

```

Usage: varcall <-s|-v> <-f REF> [options] bam1 [bam2...]
Version: 0.95.794 (BETA)

Either outputs summry stats for the list of files, or performs variant calling.  
(You can specify -s and -v, in which case both are done).

Options (later options override earlier):

-s          Calculate statistics
-v|version  Calculate variants bases on supplied parameters (see -S)
-f          Reference fasta (required if using bams, ignored otherwise)
-m          Min locii depth (1)
-a          Min allele depth (2)
-p          Min allele pct by quality (0)
-q          Min qual (3)
-Q          Min mapping quality (0)
-b          Min pct balance (strand/total) (0)
-D FLOAT    Max duplicate read fraction (depth/length per position) (1)
-d FLOAT    Minimum diversity (CV from optimal depth) (0.25)
-G FLOAT    Minimum agreement (Weighted CV of positional variation) (0.25)
-0          Zero out all filters, set e-value filter to 1, report everything
-B          If running from a BAM, turn off BAQ correction (false)
-R          Homopolymer repeat indel filtering (8)
-e FLOAT    Alpha filter to use, requires -l or -S (.05)
-g FLOAT    Global minimum error rate (default: assume phred is ok)
-l INT      Number of locii in total pileup used for bonferroni (1 mil)
-x CHR:POS  Output this pos only, then quit
-S FILE     Read in statistics and params from a previous run with -s (do this!)
-A ANNOT    Calculate in-target stats using the annotation file (requires -o)
-o PREFIX   Output prefix (works with -s or -v)
-F files    List of file types to output (var, varsum, eav, vcf)

Extended Options

--pcr-annot   BED      Only include reads adhering to the expected amplicons
--stranded    TYPE     Can be FR (the default), FF, FR.  Used with pcr-annot
--diversity|d FLOAT    Alias for -d
--agreement|G FLOAT    Alias for -G
--no-indels            Ignore all indels

Input files

Files must be sorted bam files with bai index files available.
Alternatively, a single pileup file can be supplied.

Output files

Varcalls go to stdout.  Stats go to stdout, or stderr if varcalling too

If an output prefix is used, files are created as follows:
   PREFIX.var         Variant calls in tab delimited 'varcall' format
   PREFIX.eav         Variant calls in tab delimited 'ea-var' format
   PREFIX.cse         Variant calls in tab delimited 'varprowl' format
   PREFIX.vcf         Variant calls, in vcf format
   PREFIX.varsum      Summary of variant calls
   PREFIX.tgt.var     On-target version of .var
   PREFIX.tgt.cse     On-target version of .cse
   PREFIX.tgt.varsum  On-target version of .varsum

Stats Output:

Contains mean, median, quartile information for depth, base quality, read len,
mapping quality, indel levels. Also estimates parameters suitable for
variant calls, and can be passed directly to this program for variant calls

If an output prefix is used, files are created as follows:

   PREFIX.stats       Stats output
   PREFIX.noise       Non-reference, non-homozygous allele summary
   PREFIX.xnoise      Like noise, but with context-specific rates


```