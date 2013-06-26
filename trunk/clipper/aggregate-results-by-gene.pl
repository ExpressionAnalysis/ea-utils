#!/usr/bin/perl

#Copyright (c) 2012 Erik Aronesty (erik@q32.com)
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in
#all copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#THE SOFTWARE.
#
#ALSO, IT WOULD BE NICE IF YOU LET ME KNOW YOU USED IT.

use strict;
use Getopt::Long;
use EA;
use Data::Dumper;

my %opt;
GetOptions(\%opt, "gtf=s") || die usage();
die usage() if !$opt{gtf};

# open gtf
open(G, $opt{gtf}) || die("$opt{gtf}: $!\n");
my $f  = shift @ARGV;

# open counts
open(IN, $f) || die "$f: $!\n";

# read gtf into hash
my %gs;
my %tids;
while(<G>) {
    my @fds = split /\t/, $_;
    my ($g) = $fds[8]=~/gene_id "([^"]+)"/;
    my ($t) = $fds[8]=~/transcript_id "([^"]+)"/;
    $gs{$t}=$g;
    $tids{$g}{$t}=1;
}

die "gtf file has no genes\n" unless %gs;
my $head1 = <IN>; chomp $head1;

my $format = $head1 =~ /est_counts/                 ? "express" : 
             $head1 =~ /^\S*:\S+\t\d+\t\d+\t\d/     ? "coverage-matrix" : 
             undef;

die "File is not known express 1.x, or sam-stats coverage-matrix output\n" unless $format;

close IN;

for (keys(%tids)) {
    $tids{$_}=join ',', keys(%{$tids{$_}});
}

my @fields = ();
if ($format eq 'coverage-matrix') {
    @fields = (fields=>["transcript_id", "length", "count", "coverage", "skew"]);
}

my %agg;
my $tot_count;
xsvparse($f, sub => sub {
    my ($d) = @_;
    if ($format eq 'express') {
        my $g = $gs{$d->{target_id}};
        return unless $g;
        $agg{$g}{fpkm}+=$d->{fpkm};
        $agg{$g}{eff_length}+=$d->{eff_length}*($d->{est_counts}+1);
        $agg{$g}{length}+=$d->{length}*($d->{est_counts}+1);
        $agg{$g}{est_counts}+=$d->{est_counts};
        $agg{$g}{eff_counts}+=$d->{eff_counts};
        $agg{$g}{junk_counts}+=$d->{est_counts}+1;
        $agg{$g}{fpkm_conf_high}+=$d->{fpkm_conf_high};
        $agg{$g}{fpkm_conf_low}+=$d->{fpkm_conf_low};
        $agg{$g}{transcript_ids}=$tids{$g};
    } else {
        $tot_count+=$d->{count};
        my $g = $gs{$d->{transcript_id}};
        return unless $g;
        $agg{$g}{fpkm}+=1000*$d->{count}/$d->{length};
        $agg{$g}{count}+=$d->{count};
        $agg{$g}{length}+=$d->{length}*($d->{count}+1);
        $agg{$g}{transcript_ids}=$tids{$g};
        $agg{$g}{junk_counts}+=$d->{count}+1;
    }
},nocase=>1, @fields);

my @fds = qw(gene_id est_counts eff_counts length eff_length fpkm fpkm_conf_low fpkm_conf_high transcript_ids);

if ($format eq 'coverage-matrix') {
    @fds = qw(gene_id count length fpkm transcript_ids);
}
print join("\t", @fds), "\n";

for my $g (sort {$agg{$b}->{fpkm} <=> $agg{$a}->{fpkm}} keys(%agg) ) {
    $agg{$g}{gene_id}=$g;
    $agg{$g}{eff_length}/=$agg{$g}{junk_counts};
    $agg{$g}{length}/=$agg{$g}{junk_counts};
    if ($format eq 'coverage-matrix') {
        $agg{$g}{fpkm}/=($tot_count/1000000);
    }
    print join("\t", map {$agg{$g}{$_}} @fds), "\n";
}
 

sub usage {
<<EOF
usage: $0 -g gtf_file IN1 IN2 ...

Aggregates isoform counts by the gene information in a GTF file gene info

EOF
}

