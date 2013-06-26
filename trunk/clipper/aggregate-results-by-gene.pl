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
use IO::File;
use Data::Dumper;
use Carp qw(croak);

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

Works with express-1.X output, or ea-utils/sam-stats coverage-matrix output.

EOF
}


# this ia a really cool function that turns a csv file into a hash table
# it's not "complete" but it works well... and at least i understand it
# usage: xsvparse($file, option1=>val, option1=>val)
#
# all options are optional... by defualt it guesses the quote char, and the delimiter, and returns an array of hashes
# if you specify a key, then it will return a hash of hashes indexed by that key
# for example xsvparse("Summary.tsv", key=lane);

# OPTIONS:
#   nocase  BOOL      : lowercases all keys and field names
#   skip    INT|REGEX : skips lines of file before reading
#   fields  AREF      : reference to an array of field names, for files without headers
#   key     NAME      : if specified, you get a hash of hashes indexed by a key field
#   delim   TEXT      : if specified, the delimiter to use (otherwise autodetect)
#   quot    TEXT      : if specified, the quot char to use (otherwise autodetect)
#   multi   BOOL      : if specified, each keyed entry is an arrayref of all matches
#   sub     FUNCREF   : if specified, each row is passed to the funcref


sub xsvparse {
    my ($file, %op) = @_;

    my $in = new IO::File;
    my $fin = $file =~ /\.gz$/ ? "gunzip -c $file|" : $file;
    open ($in, $fin) || return undef;

    lcasehash(\%op);

    my $l1;

    if (!$op{delim}) {
        $l1 = <$in>;
        $l1 =~ s/[\r\n]+$//o;
    }

    my $skip_me = sub {
        if ($op{skip}) {
            if (ref($op{skip}) eq 'Regexp') {
                do {
                    $l1 = <$in>; $l1 =~ s/[\r\n]+$//o;
                } while ($l1 =~ $op{skip});
            } else {
                my $skip = $op{skip};
                while ($skip > 0) {
                    $l1 = <$in>; $l1 =~ s/[\r\n]+$//o;
                    --$skip;
                }
            }
        }
    };

    &$skip_me();

    return undef unless $l1 || $op{delim};

    local $_;

    if (!$op{delim}) {
        my $n = 0;
        my $d;
        for ("\t", ",", "|") {
            my $c = $_ eq "\t" ? $l1 =~ tr/\t// : $_ eq "," ?  $l1 =~ tr/,// : $l1 =~ tr/|//;
            if ($c > $n) {$d = $_; $n = $c};
        }
        $op{delim} = $d if $n > 0;
        croak "xsvparse can't determine delimiter for $file" if $n == 0;
        if ($op{fields}) {
            close($in);
            open($in, $fin);
            &$skip_me();
        }
    }

    my $rxm = undef;
    if (!$op{quot}) {
        my $n = $l1 =~ tr/\"//;
        if ($n > 1) {
            $op{quot} = '"';
        }
    }
    if ($op{quot}) {
        # quoted csv parser... slower, but works fine
        $rxm = qr{$op{quot}([^$op{quot}\\]*(?:\\.[^$op{quot}\\]*)*)$op{quot}$op{delim}?|([^$op{delim}]+)$op{delim}?|$op{delim}};
    }
    if ($op{fields}) {
        if (!ref($op{fields}) eq 'ARRAY') {
            die "list of field names should be an array";
        }
    } else {
        if (!$op{fields}) {
            $l1 = <$in> if ! defined $l1;
            if ($rxm) {
                @{$op{fields}} = ();
                push(@{$op{fields}} ,$+) while $l1 =~ m/$rxm/gx;
            } else {
                $op{fields} = [split($op{delim}, $l1)];
                map {s/\s+$//} @{$op{fields}};
            }
        }
        if ($op{rnames}) {
            $op{fields} = ['<key>', @{$op{fields}}] unless $op{fields}->[0] eq '';
            $op{fields}->[0]=$op{key}='<key>';
        }
        if ($op{nocase}) {
            for (@{$op{fields}}) {
                $_=lc($_);
            }
        }
    }
    if ($op{nocase} && $op{key}) {
        $op{key} = lc($op{key});
    }

    no warnings;
    my $r;
    my @d;
    my @d;
    while (<$in>) {
        s/[\r\n]+$//o;
        if ($rxm) {
            @d = ();
            push(@d ,$+) while m/$rxm/gx;
            push(@d,undef) if substr($_, -1,1) eq $op{delim};
        } else {
            @d = split($op{delim}, $_);
        }
        my %h;
        for (my $i = 0; $i < @{$op{fields}}; ++$i) {
            $h{$op{fields}->[$i]}=$d[$i];
        }
        if ($op{sub}) {
            last if &{$op{sub}}(\%h) eq 'last';
        } elsif ($op{key}) {
            if ($op{nocase}) {
                if ($op{multi}) {
                    push @{$r->{lc($h{$op{key}})}}, \%h;
                } else {
                    $r->{lc($h{$op{key}})}=\%h;
                }
            } else {
                if ($op{multi}) {
                    push @{$r->{$h{$op{key}}}}, \%h;
                } else {
                    $r->{$h{$op{key}}}=\%h;
                }
            }
        } else {
            push @{$r}, \%h;
        }
    }

    return $r;
}

sub lcasehash {
    my $h = shift;
        for (keys(%$h)) {
                if (! $_ eq lc($_) ) {
                        $h->{lc($_)}=$h->{$_};
                        delete $h->{$_};
                }
        }
}

