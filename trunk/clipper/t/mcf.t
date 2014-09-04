use Test::Builder;
use Test::More;
use File::Basename qw(dirname);

require (dirname(__FILE__) . "/test-prep.pl");

$prog="$BINDIR/fastq-mcf";

@check = (
    {param=>"-l 15 $INDIR/test.fa $INDIR/test1.fq > %o:$TMPDIR/test1.out 2> %o:$TMPDIR/test1.err"},
    {param=>"-l 15 $INDIR/test.fa $INDIR/test2.fq > %o:$TMPDIR/test2.out 2> %o:$TMPDIR/test2.err", bad=>1},
    {param=>"-l 15 $INDIR/test.fa $INDIR/test1.fq -o %o:$TMPDIR/test3.out > %o:$TMPDIR/test3.err 2>&1"},
    {param=>"-l 15 -L72 -f $INDIR/test.fa $INDIR/test4.fq1 $INDIR/test4.fq2 -o %o:$TMPDIR/test4.out1 -o %o:$TMPDIR/test4.out2 > %o:$TMPDIR/test4.err 2>&1"},
    {param=>"-l 15 -H $INDIR/test.fa $INDIR/test1.fq -S -o %o:$TMPDIR/test5.out > %o:$TMPDIR/test5.err 2>&1"},
    {param=>"-l 15 $INDIR/test.fa $INDIR/test1.fq -o %o:$TMPDIR/test6.out.gz > %o:$TMPDIR/test6.err 2>&1"},
    {param=>"-0 -D 20 n/a $INDIR/test-mcf-dup.fq -o %o:$TMPDIR/test7.out > %o:$TMPDIR/test7.err 2>&1"},
    {param=>"n/a $INDIR/count.fq > %o:$TMPDIR/test8.out 2> %o:$TMPDIR/test8.err 2>&1"},
    {param=>"$INDIR/adap.fa $INDIR/test5.fq > %o:$TMPDIR/test9.out 2> %o:$TMPDIR/test9.err"},
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
