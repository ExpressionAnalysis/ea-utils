#!/usr/bin/perl

use strict;
use Data::Dumper;

while(<>) {
    next if /^#/;
    last;
}

chomp;
my @h = split /","/;

while(<>) {
    chomp;
    my @f = split /","/;
    my $i = 0;
    my %d;
    for (@h) {
        s/"//g;
        $d{lc($_)}=$f[$i++];
    }
    for (split m{///}, $d{alignments}) {
        next if /---/;
        my ($al) = split m{//};
        $al=~s/^\s+//;
        $al=~s/\s+$//;
        $al=~s/\t/ /g;
        my ($chr, $st, $en, $strand) = $al =~ /([^:]+):(\d+)-(\d+) \((.)\)/;
        --$st;
        print "$chr\t$st\t$en\t{$al}\t0\t$strand\n";
    }
}

