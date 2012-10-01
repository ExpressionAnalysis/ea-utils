#!/usr/bin/perl

# Copyright (c) 2011 Erik Aronesty (erik@q32.com)
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
# 
# ALSO, IT WOULD BE NICE IF YOU LET ME KNOW YOU USED IT.

use strict;

use Data::Dumper;

my $in = shift @ARGV;

open (IN, ($in =~ /\.gz$/ ? "gunzip -c $in|" : $in =~ /\.zip$/ ? "unzip -p $in|" : "$in")) || die "$!";
my (%mrnas, %gids, %tids, $haveexon, $gff);
while (<IN>) {
	$gff = 2 if /^##gff-version 2/;
	$gff = 3 if /^##gff-version 3/;
	next if /^#/;

	s/\s+$//;
	# 0-chr 1-src 2-feat 3-beg 4-end 5-scor 6-dir 7-fram 8-attr
	my @f = split /\t/;

	my ($chr, $fil, $typ, $beg, $end, undef, $dir, $frame, $attr) = @f;

    # most ver 2's stick gene names in the id field
    my ($tid) = $attr =~ /\bParent="?([^";]+)"?/;
    my ($id) = $attr =~ /\bID="?([^";]+)"?/;
    my ($name) = $attr =~ /\bName="?([^";]+)"?/;

	next unless $tid;

	if ($typ eq 'exon' || $typ eq 'CDS') {
	    die "no position at $typ on line $." if ! $beg;
        for (split(/,/,$tid)) {
            push @{$tids{$_}}, \@f;
        }
        $haveexon = 1 if $typ eq 'exon';
    }
	if ($typ eq 'mRNA') {
        $mrnas{$id} = $tid;
    }
	if ($typ eq 'gene' && $name) {
        $gids{$id} = $name;
    }
}

for my $tid (keys(%tids)) {
	for (sort {$a->[3] <=> $b->[3]} (@{$tids{$tid}})) {
	    my ($chr, $fil, $typ, $beg, $end, undef, $dir, $frame, $attr) = @$_;
        my ($pid) = $attr =~ /\bParent="?([^";]+)"?/;
        my ($gid) = $attr =~ /\bID="?([^";]+)"?/;
        my ($id) = $attr =~ /\bName="?([^";]+)"?/;
        my ($exn) = $gid =~ /exon:(\d+)/;
        my ($gnm) = "";

        if (!$id) {
            if ($mrnas{$pid}) {
                $id = $mrnas{$pid};
                if ($gids{$gid}) {
                    $gnm=$gids{$id};
                }
            }    
        }

        # gff3 puts :\d in exons sometimes

        $id =~ s/:\w+$//;
        $tid =~ s/:\w+$//;
        $id =~ s/-cds$//;
        $tid =~ s/^CDS_//;

        $tid = "$tid:$id" if $tid !~ /^$id/;

        my $ex = ""; 
        $ex = " exon_number \"$exn\";" if $exn > 0;

# for some this is needed... not for others?
#        --$beg; --$end;
        print "$chr\t$fil\t$typ\t$beg\t$end\t0\t$dir\t$frame\tgene_id \"$id\"; transcript_id \"$tid\";$ex\n";
        print "$chr\t$fil\texon\t$beg\t$end\t0\t$dir\t$frame\tgene_id \"$id\"; transcript_id \"$tid\";$ex\n" if !$haveexon;
	}
}

close IN;
