# Introduction #

fastq-mcf attempts to:

  * Detect & remove sequencing adapters and primers
  * Detect limited skewing at the ends of reads and clip
  * Detect poor quality at the ends of reads and clip
  * Detect Ns, and remove from ends
  * Remove reads with CASAVA 'Y' flag (purity filtering)
  * Discard sequences that are too short after all of the above
  * Keep multiple mate-reads in sync while doing all of the above

# Usage #

```

Usage: fastq-mcf [options] <adapters.fa> <reads.fq> [mates1.fq ...]
Version: 1.04.636

Detects levels of adapter presence, computes likelihoods and
locations (start, end) of the adapters.   Removes the adapter
sequences from the fastq file(s).

Stats go to stderr, unless -o is specified.

Specify -0 to turn off all default settings

If you specify multiple 'paired-end' inputs, then a -o option is
required for each.  IE: -o read1.clip.q -o read2.clip.fq

Options:
    -h       This help
    -o FIL   Output file (stats to stdout)
    -s N.N   Log scale for adapter minimum-length-match (2.2)
    -t N     % occurance threshold before adapter clipping (0.25)
    -m N     Minimum clip length, overrides scaled auto (1)
    -p N     Maximum adapter difference percentage (10)
    -l N     Minimum remaining sequence length (19)
    -L N     Maximum remaining sequence length (none)
    -D N     Remove duplicate reads : Read_1 has an identical N bases (0)
    -k N     sKew percentage-less-than causing cycle removal (2)
    -x N     'N' (Bad read) percentage causing cycle removal (20)
    -q N     quality threshold causing base removal (10)
    -w N     window-size for quality trimming (1)
    -H       remove >95% homopolymer reads (no)
    -0       Set all default parameters to zero/do nothing
    -U|u     Force disable/enable Illumina PF filtering (auto)
    -P N     Phred-scale (auto)
    -R       Dont remove Ns from the fronts/ends of reads
    -n       Dont clip, just output what would be done
    -C N     Number of reads to use for subsampling (300k)
    -S       Save all discarded reads to '.skip' files
    -d       Output lots of random debugging stuff

Quality adjustment options:
    --cycle-adjust    CYC,AMT     Adjust cycle CYC (negative = offset from end) by amount AMT
    --phred-adjust    SCORE,AMT   Adjust score SCORE by amount AMT

Filtering options*:
    --[mate-]qual-mean  NUM       Minimum mean quality score
    --[mate-]qual-gt    NUM,THR   At least NUM quals > THR
    --[mate-]max-ns     NUM       Maxmium N-calls in a read (can be a %)
    --[mate-]min-len    NUM       Minimum remaining length (same as -l)
    --hompolymer-pct    PCT       Homopolymer filter percent (95)

If mate- prefix is used, then applies to second non-barcode read only

Adapter files are 'fasta' formatted:

Specify n/a to turn off adapter clipping, and just use filters

Increasing the scale makes recognition-lengths longer, a scale
of 100 will force full-length recognition of adapters.

Adapter sequences with _5p in their label will match 'end's,
and sequences with _3p in their label will match 'start's,
otherwise the 'end' is auto-determined.

Skew is when one cycle is poor, 'skewed' toward a particular base.
If any nucleotide is less than the skew percentage, then the
whole cycle is removed.  Disable for methyl-seq, etc.

Set the skew (-k) or N-pct (-x) to 0 to turn it off (should be done
for miRNA, amplicon and other low-complexity situations!)

Duplicate read filtering is appropriate for assembly tasks, and
never when read length < expected coverage.  -D 50 will use
4.5GB RAM on 100m DNA reads - be careful. Great for RNA assembly.

*Quality filters are evaluated after clipping/trimming
```

## Notes ##

Adapter file format is fasta.  You can set it to /dev/null, and pass "-f" to do skew detection only.

## Todo ##

  * When discarding one read for being "too short", it has to discard both pairs.   For a sequencing run of normal quality this is not an issue.   It should, though, write "un-mated" reads (whose mate was skipped) to a separate file.  Typically, since these read mates were poor quality, it's not really useful... but it can be for diagnostics.   I've seen runs where these provide valuable data.

  * Like any tool that does many things, fastq-mcf can be limited in it's ability to be flexible.   The biggest missing feature is for it to be able to read files that are formatted like it's stderr output, and use them to guide the process.   Given that feature, fastq-mcf would be complete.

## Notes ##

  * Default settings are probably too conservative when it comes to trimming poor quality/detecting base-skew.

  * It won't trim the "insides" of a paired-end read.   It also will no longer attempt to quality filter a barcode read.   No override for these, but I can't think of a reason to.

  * The -x percentage can be a confusing parameter.   it causes the **entire** cycle to be removed from **all** reads.... even ones without N's... if more than, say, 20% of that cycle is N's.

## Cleaning multiple files ##

Using process substitution, or named pipes, you can clean multiple fastq's in one pass.  This is useful for combining multiple MiSeq runs, or multiple lanes for example:

```
fastq-mcf \
  -o cleaned.R1.fq.gz \
  -o cleaned.R2.fq.gz \
  adapters.fa \
  <(gunzip -c uncleaned.lane1.R1.fq.gz uncleaned.lane2.R1.fq.gz;) \
  <(gunzip -c uncleaned.lane1.R2.fq.gz uncleaned.lane2.R2.fq.gz;)
```

(Many bioinformatic tools are not "stream friendly",  and some may require the "buffer" command to work.   But fastq-mcf does its own buffering internally.)