use strict;
use Test::More;

# core modules only!
use File::Basename qw(dirname basename);
use File::Compare;
use File::Path;
use File::Spec;
use File::Temp;

use strict;

our $COPY_OK=0;

use Getopt::Long qw(:config pass_through no_ignore_case);

GetOptions("X"=>\$COPY_OK);

my $testdir = dirname(__FILE__);

chdir($testdir) || die("$testdir : $!\n");

my $tempbase = "tmp";

mkdir($tempbase);

my $template = basename($0) . ".XXXXX";

# exported 
our $BINDIR = "..";
our $TMPDIR=File::Temp::tempdir($template, CLEANUP=>0, DIR=>$tempbase);
our $INDIR="in/" . basename($0); $INDIR =~ s/\.t$//;
our $OUTDIR="out/" . basename($0); $OUTDIR =~ s/\.t$//;

sub check_output {
    my ($arr) = @_;
    for my $f (@$arr) {
        my $o = File::Spec->catfile($OUTDIR,basename($f));
        my $c1 = ($f =~ /.gz$/ ? "gunzip -c '$f'" : "cat '$f'") . "|perl -pe 's/\Q$TMPDIR\E/#TMPDIR#/g' |";
        if ($COPY_OK) {
            my $cp = "$c1 " . ( $f =~ /.gz$/ ? "gzip -c" : "cat" ) . " > '$o'";
            system($cp);
        }
        my ($i1, $i2);
        my $c2 = ($o =~ /.gz$/ ? "gunzip -c '$o'" : "cat '$o'") . "|perl -pe 's/\Q$TMPDIR\E/#TMPDIR#/g' |";
        open $i1, $c1;
        open $i2, $c2;
        ok(compare($i1, $i2) == 0, "Files equal: $f == $o");
    }
}

sub run {
    my ($cmd) = @_;
    my @o = $cmd =~ m/[#%]o:(\S+)/g;

    chomp $cmd;

    $cmd =~ s/%o:(\S+)/$1/g;
       
#    warn "# $cmd\n";

    my $exit = system($cmd);

    if ($exit >> 8) {
        $exit = $exit >> 8;
    }

    return $exit, $cmd, [@o];
}

END {
    if (Test::More->builder->is_passing() && !$?) {
        warn "# removing $TMPDIR\n";
        File::Path::rmtree($TMPDIR);
    }
}


1;
