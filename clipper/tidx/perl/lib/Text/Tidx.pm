package Text::Tidx;

use 5.010001;
use strict;
use warnings;
use Carp;

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use Text::Tidx ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw(
	
) ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
	
);

our $VERSION = '0.94';

require XSLoader;
XSLoader::load('Text::Tidx', $VERSION);

# Preloaded methods go here.

sub lookup {
    if (@_ == 3) {
        return lookup_c(@_, "^");
    }
    if (@_ == 4) {
        return lookup_cr(@_, "^");
    }
}

sub query {
    my $r;
    if (@_ == 3) {
        $r = lookup_c(@_, "^");
    }
    if (@_ == 4) {
        $r = lookup_cr(@_, "^");
    }
    return () if (!$r);
    $r =~ s/^\^//;
    return split /\^/, $r;
}

sub build {
    my ($file, %op) = @_;
    croak "usage: build(file, options)\n" unless $file;
    $op{skip} = '#' if !defined($op{skip});
    $op{skip_c} = $op{skip} !~ /^\d+$/ ? $op{skip} : '';
    $op{skip_i} = $op{skip} =~ /^\d+$/ ? $op{skip} : 0;
    $op{sub_e} = $op{sub} || $op{sub_e} ? 1 : 0;
    $op{sep} = "\t" if !$op{sep};
    $op{chr} = 1 if !defined($op{chr});
    $op{beg} = 2 if !$op{beg};
    $op{end} = 3 if !$op{end};
    # one based index, consistent with command-line version
    --$op{chr};
    --$op{beg};
    --$op{end};
    # todo... for kicks: allow indexing on text only, no positions
    tidx_build($file, $op{sep}, $op{chr}, $op{beg}, $op{end}, $op{skip_i}, $op{skip_c}, $op{sub_e});
}

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

Text::Tidx - Index a delimited text file containing start-stop positions

=head1 SYNOPSIS

  use Text::Tidx;
  Text::Tidx::build("annot.txt");
  $idx = Text::Tidx->new("annot.txt");
  print $idx->query("chr1",240034);

=head1 FUNCTION

=head2 new(FILE)

Loads an index from a file.

=head2 query(CHR, POS [, END])

Query a loaded index, returning an array of text lines corresponding to the specified
chr string and integer pos.   If an end is specified, then all overlapping regions
are returned.

=head2 build(FILE [, option1=>value1, ...])

Builds an index.  Default is to index on the first 3 columns.

The following options may be used:

=over 4

=item sep

Field separator, default to a tab

=item chr

1-based index of the string key field, can be -1 for "Not applicable", default is 1

=item beg

1-based index of the field containing the start of the integer numeric range, default is 2

=item end

1-based index of the field containing the end of the integer numeric range, default is 3

=item skip

If an integer, then it is the number of rows to skip.  If it's 
a character, then skips all rows beginning with that character.
Default is '#', skipping comment chars (compatible with gffs, vcfs, etc.)

=item sub_e

If nonzero, then the "end" of the range is not included in the range, ie: 
one is subtracted from the end positions.

=back

=head1 DESCRIPTION

Text:Tidx allows you to index any text file using a key field
and range coordinates, and, later, use that index for O(log(n))
range-lookups into the file.

This was written because it was, for me significantly faster, for very large
files (>100k rows) and many searches ( > 10k), then entering
all of the information into a database and doing range querys,
even faster than SQLITE's rtree extension, or the "tabix" program
both of which are do similar things and do them rather well.

Although it was designed for chromosome, stop, start indexing,
it is not genome specific, and can index any delimited text file.

Indexes are loaded into RAM.  If you only have a few lookups
to do perl instance, this is expensive, and a database will be faster.

=head1 AUTHOR

Erik Aronesty, E<lt>earonesty@cpan.orgE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2012 by Erik Aronesty

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.1 or,
at your option, any later version of Perl 5 you may have available.


=cut
