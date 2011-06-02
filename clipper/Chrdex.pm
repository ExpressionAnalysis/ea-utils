package Chrdex;

# Licensed via the "Artistic License" 
# See: http://dev.perl.org/licenses/artistic.html
# Copyright 2011, Expression Analysis
# Author: Erik Aronesty <earonesty@expressionanalysis.com>
# "Let me know if it's useful"
#
# EXMAPLE:
# $x = Chrdex->new("CCDS_Exome_annot.txt", chr=>2, beg=>5, end=>6, skip=>1);
# $x->search(1, 153432255);
# TODO:
# work with ranges...should be easy

use Inline 'C';

use strict;
use warnings::register;

use Storable qw(store retrieve);
use Data::Dumper;
use locale; ##added by vjw to control case of reference bases

our $VERSION = '1.2.15';		# major = object interface is different, minor = new features, release=fixes/performance improve
my $FILEVER = 4;			# only increment this if existing files won't work with the new version

my $tmpb;

sub new {
	my ($class, $path, %opts) = @_;

	if (ref($class)) {
		$class = ref($class);
	}

	$opts{delim} = "\t" if ! $opts{delim};
	$opts{skip} = 0 if ! $opts{skip};
	$opts{chr} = 0 if ! $opts{chr};
	$opts{beg} = 1 if ! $opts{beg};
	$opts{end} = 2 if ! $opts{end};
	$opts{ver} = $FILEVER;

	if (! -s $path) {		# be a little more careful about this one
		if (! -s $path || -d $path) {
			die "Can't open $path.\n";
		}
	}

	# location of data store
	my $annob = $path;
	$annob =~ s/([^\/]+)$/\.$1/;
	$annob = "$annob.chrdex";

	$annob = $opts{index_path} if $opts{index_path};

	my $ref;
        my $mt = (stat($path))[9];
	# if index is new
	if (!$opts{force} && (stat($annob))[9] > $mt) {
		$ref = eval {retrieve $annob};

		# if arguments were different, then clear ref
		for (qw(delim skip chr beg end ver byref)) {
			last if !$ref;
			$ref = undef if !($ref->{_opts}->{$_} eq $opts{$_});
		}

		if ($ref) {
			# if begin != end, then type is range
			if ($ref->{_opts}->{beg} != $ref->{_opts}->{end}) {
				eval{chrdex_check($ref)};
				if ($@) {
					$ref = undef;
				}
			}
		}
	}

	if (!$ref) {
		my %locs;
		$tmpb = "$annob.$$";
		open( IN, $path ) or die "Can't open $path: $!\n";
		my $skip = $opts{skip};
		while ($skip > 0) { scalar <IN>; --$skip };

		my $pos = 0;
		$pos = tell IN if $opts{byref};
		while(<IN>) {
			my ($chr, $beg, $end);
			$_ =~ s/\s+$//;
			my @data = split /\t/;
			$chr = $data[$opts{chr}];
			$beg = $data[$opts{beg}];
			$end = $data[$opts{end}];	
			if (!(($beg+0) eq $beg)) {
				die "Invalid data in $path at line $., expected a number, got '$beg'\n";
			}
			$chr=~s/^chr//i;
			$_ = $pos if $opts{byref};

			# here's where you put the annotation info
			if ($opts{beg} == $opts{end}) {
				if ($locs{"$chr:$beg"}) {
 					push @{$locs{"$chr:$beg"}}, $_;;
				} else {
					$locs{"$chr:$beg"} = [$_];
				}
			} else {
				push @{$locs{$chr}}, [$beg+0, $end+0, [$_]];
			}
			$pos = tell IN if $opts{byref};
		}
		close IN;

		if ($opts{beg} == $opts{end}) {
			goto DONE;
		}

		# sort & cache annotation, deal with overlaps nicely
		my $i;
		for my $chr (keys(%locs)) {
			my $arr = $locs{$chr};
			@{$locs{$chr}} = sort {$a->[0]-$b->[0]} @{$locs{$chr}};
			for ($i=0;$i<$#{$arr};++$i) {
				next unless $arr->[$i+1]->[0];				# empty? skip
				if ($arr->[$i]->[1] >= $arr->[$i+1]->[0]) {		# if i overlap the next one
					# warn 1, Dumper($arr->[$i], $arr->[$i+1], $arr->[$i+2]);
				
					# frag after next	
					my $new_st = $arr->[$i+1]->[1]+1;
					my $new_en = $arr->[$i]->[1];
					my $new_ro = $arr->[$i]->[2];
			
					# TODO: store as array... string folding will save lots of space when there are many overlaps
					# but hasn't been a problem so far
					if ($arr->[$i]->[1] < $arr->[$i+1]->[1]) {
						# overlap next
						$new_st = $arr->[$i]->[1] + 1;
						$new_en = $arr->[$i+1]->[1];
						$new_ro = [@{$arr->[$i+1]->[2]}];
						$arr->[$i+1]->[1] = $arr->[$i]->[1];
						push @{$arr->[$i+1]->[2]}, @{$arr->[$i]->[2]};
					} else {
						push @{$arr->[$i+1]->[2]}, @{$arr->[$i]->[2]};
					}

					# shorten my end to less than the next's start
					$arr->[$i]->[1] = $arr->[$i+1]->[0]-1;

					if ($new_en >= $new_st) {
						# warn "NEW: $new_st $new_en $new_ro\n";

						# put the fragment where it belongs
						my $j=$i+2;
						while ($j<=$#{$arr} & $new_st > $arr->[$j]->[0]) {
							++$j;
						}
						splice(@{$arr}, $j, 0, [$new_st, $new_en, $new_ro]);
					}
					
					if ($arr->[$i]->[1] < $arr->[$i]->[0]) {
						splice(@{$arr}, $i, 1);
						--$i;
					}
					# warn 2, Dumper($arr->[$i], $arr->[$i+1], $arr->[$i+2]);
				}
			}
		}
		DONE:
		$locs{_opts} = \%opts;
		$ref = \%locs;
		store \%locs, "$tmpb";
		rename "$tmpb", "$annob";
	}

	# stuff these into the top level, since tests showed it was significantly more expensive to doubly-reference

	$ref->{_type}='C';
	$ref->{_type}='I' if ($opts{beg} == $opts{end});

	if (($ref->{_type} eq 'C')) {
		chrdex_check($ref);
	}

	if ($opts{byref}) {
		# only storing pointers to file, not whole record
		require IO::File;
		$ref->{_byref}=1;
		$ref->{_refh} = new IO::File;
		open($ref->{_refh}, $path) || die "Can't open $path\n";
	}

	bless $ref, $class;

	return $ref;
}

END {
	unlink("$tmpb.chrdex");
}

sub search {
	return join "\n", query(@_);
}

sub query {
	my ($self, $chr, $loc, $loc2) = @_;
	$chr=~s/^chr//io;
	my $type = $self->{_type};
	my $list;
	if ($type eq 'C') {
		if ($loc2) {
			$list = chrdex_search_range($self, $chr, $loc, $loc2);
			return () if ! defined $list;
			my $prev;
			for (my $i = 0; $i < @$list; ++$i) {
				if ($prev eq $list->[$i]) {
					splice @$list, $i, 1;
					--$i;
				}
				$prev = $list->[$i];
			}
			return @$list if !($self->{_byref});
		} else {
			$list = chrdex_search($self, $chr, $loc);	
			return () if ! defined $list;
			return @$list if !($self->{_byref});
		}
	} elsif ($type eq 'I') {
		if ($loc2) {
			# if only this wasn't a hash.... sheesh
			my %hv;
			for (my $i=$loc;$i<$loc2;++$i) {
				$list = $self->{"$chr:$loc"};
				next unless defined $list;
				for (@$list) {
					$hv{$_}=1;
				}
			}
			$list = [keys(%hv)];
		} else {
			$list = $self->{"$chr:$loc"};
			return () if ! defined $list;
		}
	}

	if ($self->{_byref}) {
		for (@$list) {
			seek $self->{_refh}, $_, 0;
			$_=$self->{_refh}->getline();
			chomp $_;
		}
	}
	return @$list;
}

1;

__DATA__
__C__

bool get_sten(AV *arr, int i, int *st, int*en);
SV * av_fetch_2(AV *arr, int i, int j);
int chrdex_search_n(AV *arr, SV *schr, SV* sloc);

void chrdex_search(SV *self, SV *schr, SV* sloc) {
	SV *roi;
        AV *arr;
        SV **pav;
        HV *map= (HV*) SvRV(self);
        char *chr = SvPV_nolen(schr);
        int loc = SvIV(sloc);

        pav = hv_fetch(map, chr, strlen(chr), 0);

        if (!pav)
                return;

        arr = (AV*) SvRV(*pav);

	int i = chrdex_search_n(arr, schr, sloc);
	if (i >=0) {
                roi = av_fetch_2(arr, i, 2);
                if (!roi)
                        return;

                Inline_Stack_Vars;
                Inline_Stack_Reset;
                Inline_Stack_Push(sv_2mortal(newSVsv(roi)));
                Inline_Stack_Done;
                Inline_Stack_Return(1);
	}
	return;
}

void chrdex_search_range(SV *self, SV *schr, SV* sloc, SV* eloc) {
        SV *roi = NULL;
        AV *arr;
        SV **pav;
        HV *map= (HV*) SvRV(self);
        char *chr = SvPV_nolen(schr);
        int isloc = SvIV(sloc);
        int ieloc = SvIV(eloc);

        pav = hv_fetch(map, chr, strlen(chr), 0);

        if (!pav)
                return;

        arr = (AV*) SvRV(*pav);

        int i = chrdex_search_n(arr, schr, sloc);
        int j = chrdex_search_n(arr, schr, eloc);

	if (!i) i=j;
	if (!j) j=i;
        if (i >=0) {
		int x;
		char *rx;
		AV *rav=NULL;
		for (x=i; x<=j; ++x) {
			int st, en;
			if (get_sten(arr, x, &st, &en)) {
				if (ieloc >= st && isloc <= en) {
					SV* ret = av_fetch_2(arr, x, 2);
					if (ret) {
						int z;
						if (!rav) rav = newAV();
						for (z=0;z<=av_len(SvRV(ret));++z) {
							SV ** s = av_fetch(SvRV(ret), z, 0);
							if (s) {
								SvREFCNT_inc(*s);
								av_push(rav, *s);
							}
						}
					}
				}
			}
		}

                if (!rav)
                        return;

                Inline_Stack_Vars;
                Inline_Stack_Reset;
                Inline_Stack_Push(newRV_noinc((SV*)rav));
                Inline_Stack_Push(roi);
                Inline_Stack_Done;
                Inline_Stack_Return(1);
        }
        return;
}

int chrdex_search_n(AV *arr, SV *schr, SV* sloc) {
	int b=0, t, i, st, en;
	char *chr = SvPV_nolen(schr);
	int loc = SvIV(sloc);

	b = 0;
	t = av_len(arr);

	if (t <= b) {
		get_sten(arr, i=0, &st, &en);
	} else {
            while (t > b) {
                i = (t+b)/2;
		if (!get_sten(arr, i, &st, &en))
			return;

                if ((i == b) || (i == t)) 
			break;
                if (loc > en) {
                        b = i;
                } else if (loc < st) {
                        t = i;
                } else {
                        break;
                }
            }
	}

//	printf("chr:%s loc: %d, st: %d en: %d i: %d t: %d b: %d\n", chr, loc, st, en, i, t, b);

	if (loc < st) {
		--i;
		if (i < 0 || !get_sten(arr, i, &st, &en))
			return;
	} else if (loc > en) {
		++i;
		if (!get_sten(arr, i, &st, &en))
			return;
	}

        if (loc >= st && loc <= en) {
		return i;
        }

        return -1;
}

// doubly indexed array ... fetch 1, fetch 2
SV * av_fetch_2(AV *arr, int i, int j) {
        SV **pav;

        if (!(pav = av_fetch(arr,i,0)))
                return &PL_sv_undef;

        arr = (AV*) SvRV(*pav);

        if (!(pav = av_fetch(arr,j,0)))
                return &PL_sv_undef;

        return *pav;
}

bool get_sten(AV *arr, int i, int *st, int*en) {
	SV **pav;

	if (!(pav = av_fetch(arr,i,0)))
		return 0;

	arr = (AV*) SvRV(*pav);

	if (!(pav = av_fetch(arr,0,0)))
		return 0;
	*st = SvIV(*pav);

	if (!(pav = av_fetch(arr,1,0)))
		return 0;

	*en = SvIV(*pav);

	return 1;
}

void chrdex_check(SV *annoR) {
	HE *he;		// hash entry
	SV *ent;	// hash value
	HV *hv;		// annotation hash table
	AV *av;		// array in hash
	SV **v;		// entry in array
	char *key;
	int len;

        if (!SvRV(annoR))
                croak("annotation hash must be a reference");

        annoR = SvRV(annoR);
        if (SvTYPE(annoR) != SVt_PVHV)
                croak("annotation array must be a hash ref");

	hv = (HV *) annoR;
	if (!hv_iterinit(hv)) { 
                croak("empty hash, fix that in perl");
	}

	he = hv_iternext(hv);
	ent = hv_iternextsv(hv, &key, &len);
	while (key && *key == '_') {
		ent = hv_iternextsv(hv, &key, &len);
	}
	if ( SvTYPE(ent) != SVt_RV || (SvTYPE(SvRV(ent)) != SVt_PVAV) ) {
		croak("each entry in the annotation hash must be a reference to an array");
	}

	av = (AV*) SvRV(ent);
	v = av_fetch(av, 0, 0);
	if (!v) {
		croak("no empty annotation arrays, please");
	}
	
	if ( SvTYPE(*v) != SVt_RV || (SvTYPE(SvRV(*v)) != SVt_PVAV) ) {
		croak("each entry in the array should contain a start and end region");
	}

	// ok.... reference should be safe enough not to segfault later
}


