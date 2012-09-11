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

our $VERSION = '0.91';

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
    $op{chr} = 0 if !defined($op{chr});
    $op{beg} = 1 if !$op{beg};
    $op{end} = 2 if !$op{end};
    tidx_build($file, $op{sep}, $op{chr}, $op{beg}, $op{end}, $op{skip_i}, $op{skip_c}, $op{sub_e});
}

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

Text::Tidx - Perl extension for blah blah blah

=head1 SYNOPSIS

  use Text::Tidx;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for Text::Tidx, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

None by default.



=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

A. U. Thor, E<lt>earonesty@E<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2012 by A. U. Thor

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.1 or,
at your option, any later version of Perl 5 you may have available.


=cut
