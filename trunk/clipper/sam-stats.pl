#!/usr/bin/perl

# Copyright (c) 2011 Expression Analysis / Erik Aronesty
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

use strict;
use Getopt::Long;
use Data::Dumper;


my ($cq);

# summary stats
my $n = 0;		# numreads
my (%mapb, $mapb);	# mapped bases
my ($lenmin, $lenmax, $lensum, $lenssq, $mapsum, $mapssq, $nmlen, $nmc, $insum, $dlsum, $nmsum, $nmnz);
my ($nbase, $qualmax, $qualmin, $qualsum, $qualssq, @basecnt);
my ($rev, $for);
my ($mapc, @mapq, @mat);

# per-read
my ($id, $bits, $nmo, $mapq, $cig, $l, $q, $r, $mat, $nm, @fdx);

my %opt;
GetOptions(\%opt, "help", "bases=i", "make|M", "count=i", "ext|x=s") || die usage();
if ($opt{help}) { print usage(); exit 1;}

my $multi = 1 if @ARGV > 1;
$opt{ext}='stats' if $multi && ! $opt{ext};
$multi = 1 if $opt{ext};

$qualmin = 100000;
$lenmin = 100000;

my $QSAMP=$opt{bases} ? $opt{bases} : 5000;	# number of reads to track per-base quality stats (all other stats are for all reads)
my $MAXN=$opt{count} ? $opt{count} : 2000000000;	# number of reads to track per-base quality stats (all other stats are for all reads)

if (!@ARGV) {
	$ARGV[0] = '-';
}
while (@ARGV) {
my $in = shift @ARGV;
if ($in && !($in eq '-')) {
	open(IN, ($in =~ /\.gz$/ ? "/bin/gunzip -c $in|" : $in)) || die "Can't open $in: $!\n";
} else {
	die "Can't read from stdin in multiple mode\n" if ($multi);
	*IN=*STDIN;
}

*OUT=*STDOUT;
my $out;
if ($multi) {
	$out = $in;
	$out =~ s/\.gz$//;
	$out .= ".$opt{ext}";

	if ($opt{make}) {
		next if -s $out && (stat($out))[9] > (stat($in))[9];
		warn "+$in>$out\n";
	}

	open(OUT, ">$out.tmp") || die "Can't open $out\n";
}

	(%mapb, $mapb, @basecnt, $lenmin, $lenmax, $lensum, $lenssq, $mapsum, $mapssq, $nmlen, $nmc, $insum, $dlsum, $nmsum, $nmnz, $nbase, $qualmax, $qualmin, $qualsum, $qualssq, $rev, $for, $mapc, @mapq, @mat, $n) = (undef) x 100;

	%mapb=();
if ($in =~ /\.bam$/) {
	close IN;

	require Bio::DB::Sam;
	my $bam = Bio::DB::Bam->open($in);
	my $header = $bam->header;
	# reset stats
	while (my $a = $bam->read1()) {
		++$n;
		$cig = $a->cigar_str;
		next if ($cig eq '');
		$nmo=$a->tid;
		$mapq = $a->qual;
		$bits = $a->flag;
		$l = $a->query->length;
		$mat = $a->isize;
		($nm) = $a->get_tag_values("NM");
		if ($mapc < $QSAMP) {
			$r = $a->query->dna;
			$q = join '', map {chr($_+33)} unpack('C*',$a->_qscore);
		}
		stats();
		last if $n > $MAXN;
	}
} else {
	while (my $i = <IN>) {
		next if $i =~ /^\@/;
		++$n;
		($id, $bits, $nmo, undef, $mapq, $cig, undef, undef, $mat, $r, $q, @fdx) = split /\t/, $i;
		next if ($cig eq '*');
		$nm = '';
		for (@fdx) {
			if (s/^NM:i://) {
				$nm = $_;
				last;
			}
		}
		$l = length($r);
		next unless $nmo;
		stats();
		last if $n > $MAXN;
	}
}

$QSAMP=$mapc if $mapc < $QSAMP;

@mapq = sort {$a <=> $b} @mapq;
@mat = sort {$a <=> $b} @mat;

# autodetect phred
my $phred = 64;
$phred = 33 if ($qualmin < 64);

printf OUT "reads\t%d\n", $n;
printf OUT "phred\t%d\n", $phred if $nbase > 0;
printf OUT "mapped reads\t%d\n", $mapc;
printf OUT "mapped bases\t%d\n", $mapb;

if ($mapc > 0) {
	printf OUT "foward\t%d\n", $for;
	printf OUT "reverse\t%d\n", $rev;

	if ($lenmax != $lenmin) {
		printf OUT "len max\t%.4f\n", $lenmax;
		if ($n > 0 ) {
			printf OUT "len mean\t%.4f\n", $lensum/$mapc;
			if ($n > 1 ) {
				printf OUT "len stdev\t%.4f\n", stdev($mapc, $lensum, $lenssq);
			}
		}
	} else {
		printf OUT "len max\t%d\n", $lenmax;
	} 

	printf OUT "mapq mean\t%.4f\n", ($mapsum/$mapc);
	printf OUT "mapq stdev\t%.4f\n", stdev($mapc, $mapsum, $mapssq);
	printf OUT "mapq Q1\t%.4f\n", quantile(\@mapq, .25);
	printf OUT "mapq median\t%.4f\n", quantile(\@mapq, .50);
	printf OUT "mapq Q3\t%.4f\n", quantile(\@mapq, .75);

	if ($nmlen > 0) {
		printf OUT "snp rate\t%.6f\n", ($nmsum/$nmlen);
		if ($insum > 0) {
			printf OUT "ins rate\t%.6f\n", ($insum/$nmlen);
		}

		if ($dlsum > 0) {
			printf OUT "del rate\t%.6f\n", ($dlsum/$nmlen);
		}
		printf OUT "pct mismatch\t%.4f\n", 100*($nmnz/$nmc);
	}

	if (@mat > 0) {
		my $p10=quantile(\@mat,.10);
		my $p90=quantile(\@mat,.90);
		my ($matc, $matsum, $matssq);
		for (@mat) {
			if ($_ >= $p10 && $_ <= $p90) {
				++$matc;
				$matsum+=$_;
				$matssq+=$_*$_;
			}
		}
		printf OUT "insert mean\t%.4f\n", ($matsum/$matc);
		printf OUT "insert stdev\t%.4f\n", stdev($matc, $matsum, $matssq);
		printf OUT "insert Q1\t%.4f\n", quantile(\@mat, .25);
		printf OUT "insert median\t%.4f\n", quantile(\@mat, .50);
		printf OUT "insert Q3\t%.4f\n", quantile(\@mat, .75);
	}

	if ($nbase > 0) {
		printf OUT "base qual mean\t%.4f\n", ($qualsum/$nbase)-$phred;
		printf OUT "base qual stdev\t%.4f\n", stdev($nbase, $qualsum, $qualssq);
	}
}

if ($nbase) {
	printf OUT "%%A\t%.4f\n", 100*$basecnt[0]/$nbase;
	printf OUT "%%C\t%.4f\n", 100*$basecnt[1]/$nbase;
	printf OUT "%%G\t%.4f\n", 100*$basecnt[2]/$nbase;
	printf OUT "%%T\t%.4f\n", 100*$basecnt[3]/$nbase;
	if ($basecnt[4] > 0) {
		printf OUT "%%N\t%.4f\n", 100*$basecnt[4]/$nbase 
	}
}

if (%mapb > 1 && %mapb <= 1000) {
	for my $k (sort(keys(%mapb))) {
		my $v = $mapb{$k};
		my ($n) = $k =~ m/([\w.-]+)/;
		printf OUT "%%$n\t%.2f\n", 100*$v/$mapb;
	}
}

close OUT;
rename("$out.tmp", $out);
}

sub stats
{
	my $strand = 1 - 2 * ($bits & 16);

	++$mapc;

	$lenmax = $l if $l > $lenmax;
	$lenmin = $l if $l < $lenmax;

	$lensum += $l;
	$lenssq += $l*$l;

	if ($bits & 16) {
		$rev += 1;
	} else {
		$for += 1;
	}

	$mapsum += $mapq;
	$mapssq += $mapq*$mapq;

	push @mapq, $mapq+0;
	if (!($nm eq '')) {
		$nmsum += $nm;
		$nmnz += 1 if $nm > 0;
		++$nmc;
		$nmlen += $l;
		for ($cig =~ /D(\d+)/g) {
			$dlsum+=$_;
		}
		for ($cig =~ /I(\d+)/g) {
			$insum+=$_;
		}
	}

	if (%mapb <= 1000) {
		$mapb{$nmo}+=$l; 
	}
	
	$mapb+=$l;

	if ($mat > 0) {
		push @mat, $mat+0;
	}

	if ($mapc < $QSAMP) {
		# shorter length
		for (my $i = 0; $i < $l; ++$i) {
			$cq = ord(substr($q, $i, 1));
			++$nbase;
			$qualmax=$cq if $cq > $qualmax;
			$qualmin=$cq if $cq < $qualmin;
			$qualsum+=$cq;
			$qualssq+=$cq*$cq;
			my $cb = substr($r, $i, 1);
			$cb =~ tr/ACGTN/01234/;
			++$basecnt[$cb];
		}
	}
}

# copied from EA.pm

sub stdev($ $ $) {
        my ($cnt, $sum, $ssq) = @_;
        sqrt(($cnt*$ssq-($sum*$sum)) / ($cnt*($cnt-1)));
}

sub quantile {
        my ($a,$p) = @_;
        my $l = scalar(@{$a});
        my $t = ($l-1)*$p;
        my $v=$a->[int($t)];
        if ($t > int($t)) {
                return $v + $p * ($a->[int($t)+1] - $v);
        } else {
                return $v;
        }
}

sub usage {
return <<EOF
sam-stats [options] <sam or bam file> [<multiple files]

Options: 
  -h                : this help
  -b                : base sample size (5000)
  -c                : read sample count (0)
  -x		    : extension (stats)
  -M		    : make newer only

OUTPUT:

If one file is specified, then the output is to standard out.  If
multiple files are specified, or if the -x option is supplied,
the output file is <filename>.<ext>.  Default extension is 'stats'.

Complete Stats:

  <STATS>           : mean, max, stdev, median, Q1 (25 percentile), Q3
  reads             : # of entries in the file
  phred             : phred scale used
  mapped reads      : number of aligned reads
  mapped bases      : total of the lengths of the aligned reads
  forward           : number of forward-aligned reads
  reverse           : number of reverse-aligned reads
  snp rate          : mismatched bases / total bases
  ins rate          : insert bases / total bases
  del rate          : deleted bases / total bases
  pct mismatch      : percent of reads that have mismatches
  len <STATS>       : read length stats, ignored if fixed-length
  mapq <STATS>      : stats for mapping qualities
  insert <STATS>    : stats for insert sizes
  %<CHR>            : percentage of mapped bases per chromosome (use to compute coverage) 

Subsampled stats:
  base qual <STATS> : stats for base qualities
  %A,%T,%C,%G	    : base percentages
EOF
}

