use Test::Builder;
use Test::More;
use File::Basename qw(dirname);

require (dirname(__FILE__) . "/test-prep.pl");

$prog="$BINDIR/fastq-join";

@check = (
    {param=>"$INDIR/test-m1.fq $INDIR/test-m2.fq -o $TMPDIR/test-m. > %o:$TMPDIR/test-m.out 2>&1 #o:$TMPDIR/test-m.join"},
    {param=>"$INDIR/phred.1.fq $INDIR/phred.2.fq -o $TMPDIR/test-phred. > %o:$TMPDIR/test-phred.out 2>&1 #o:$TMPDIR/test-phred.join"},
    {param=>"-p 20 -m 5 $INDIR/test-ov-a.1.fq $INDIR/test-ov-a.2.fq -o $TMPDIR/test-ov-a. -x > %o:$TMPDIR/test-ov-a.out 2>&1 #o:$TMPDIR/test-ov-a.join"},
    {param=>"-p 20 -m 5 $INDIR/test-ov-b.1.fq $INDIR/test-ov-b.2.fq -o $TMPDIR/test-ov-b. -x > %o:$TMPDIR/test-ov-b.out 2>&1 #o:$TMPDIR/test-ov-b.join"},
    {param=>"-p 20 -m 5 $INDIR/test-ov-a.1.fq $INDIR/test-ov-a.2.fq -o $TMPDIR/test-nov-a. > %o:$TMPDIR/test-nov-a.out 2>&1 #o:$TMPDIR/test-nov-a.join"},
    {param=>"-p 20 -m 5 $INDIR/test-ov-b.1.fq $INDIR/test-ov-b.2.fq -o $TMPDIR/test-nov-b. > %o:$TMPDIR/test-nov-b.out 2>&1 #o:$TMPDIR/test-nov-b.join"},
);

my $id=0;
for (@check) {
    ++$id;
    my %d = %{$_};
    $cmd = "$prog $d{param}";
    my ($exit, $ncmd, $files) = run($cmd);
    if ($d{bad}) {
        ok($exit != 0, "test$id worked ($ncmd)");
    } else {
        ok($exit == 0, "test$id worked ($ncmd)");
    }

    check_output($files);
}

done_testing();
