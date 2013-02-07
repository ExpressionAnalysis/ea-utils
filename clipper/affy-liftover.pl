#!/usr/bin/perl

use strict;
use Getopt::Long;

my %opt;
die usage() unless GetOptions(\%opt, "csv=s", "chain=s", "output=s");
die usage() unless $opt{chain} && $opt{csv};

sub usage {
<<EOF

usage: $0 --csv=affy.csv --chain=lifover.chain.gz

options:

-o  PREFIX (defaults to CSVBASENAME-lift)

Output: PREFIX.csv

Also generates temp files with PREFIX.

EOF
}

my $out = $opt{csv};
$out =~ s/\.csv$//i;
$out = "$out-lift";
$out = $opt{output} if $opt{output};
$out =~ s/\.csv$//;

system("perl /opt/scripts/affy-csv-to-bed.pl $opt{csv} > $opt{csv}.bed") && die;
system("./liftOver -minMatch=0.75 $opt{csv}.bed $opt{chain} ${out}.bed ${out}.unchain") && die;

my $unmap = 0+`wc -l $out.unchain`;

print "warning\tlost $unmap annotations during liftover\n";

open (IN, "${out}.bed") || die "${out}.bed:$!";

my %map;
while (<IN>) {
# liftover coverts spaces to tabs in an unusual way.... fix this!
    chomp;

    die unless s/{([^{}]+)}/X/;

    my ($orig) = $1; 
    $orig =~ s/\t/ /g;
    
    my ($chr, $st, $en, undef, undef, $str) = split /\t/;
    die unless $str =~ /\+|-/;

    ++$st;
    my $new = "$chr:$st-$en ($str)";
    $map{$orig}=$new;    
}

open (IN, $opt{csv}) || die "$opt{csv} : $!";

while(<IN>) {
    print $_ && next if /^#/;
    last;
}
chomp;
my @h = split /","/, $_;
my $i_al=-1;
for (@h) {
    ++$i_al;
    last if /^alignments/i;
}
die "invalid CSV file, no alignments header" unless $i_al > 0;

open OUT, ">$out.csv";

while(my $l=<IN>) {
    print $l && next if /^#/;
    my @f = split /","/, $l;
    my @sub;
    for my $al (split m{///}, $f[$i_al]) {  
        my ($orig) = split m{//}, $al;
        $orig=~s/^\s+//;
        $orig=~s/\s+$//;
        $orig=~s/\t//g;
        my $new=$map{$orig};
        if ($new) {
            my $tmp=$al;    
            die unless $tmp=~ s/\Q$orig\E/$new/;
            push @sub, [$al, $tmp];
        }
    }
    for (@sub) {
        die unless $l=~s/\Q$_->[0]\E/$_->[1]/;
    }
    print OUT $l;
}

