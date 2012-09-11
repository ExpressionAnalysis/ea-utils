# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl Text-Tidx.t'

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More tests => 2;
use Cwd;

BEGIN { use_ok('Text::Tidx') };

#########################

# Insert your test code below, the Test::More module is use()ed here so read
# its man page ( perldoc Test::More ) for help writing this test script.

my ($subdir) = $0 =~ /(.*)[\/\\]/;

Text::Tidx::build("$subdir/annot.txt");

my $cytoTab=Text::Tidx->new("$subdir/annot.txt");

@res = $cytoTab->query("chr1", 2300001);
is($res[0],'chr1	2300000	5400000	p36.32	gpos25');

