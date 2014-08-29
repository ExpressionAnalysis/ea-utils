use Test::Builder;
use Test::More;
use File::Basename qw(dirname);

require (dirname(__FILE__) . "/test-prep.pl");

$prog="$BINDIR/fastq-multx";

@check = (
    {param=>"-l $INDIR/master-barcodes.txt $INDIR/mxtest_2.fastq $INDIR/mxtest_1.fastq $INDIR/mxtest_3.fastq -o n/a -o $TMPDIR/mxout_%_1.fq -o $TMPDIR/mxout_%_2.fq > %o:$TMPDIR/test1.out 2> %o:$TMPDIR/test1.err"},
    {param=>"-l $INDIR/master-barcodes.txt $INDIR/mxtest_2.fastq $INDIR/mxtest_1.fastq $INDIR/mxtest_3.fastq -o n/a -o $TMPDIR/mxout_%_1.fq.gz -o $TMPDIR/mxout_%_2.fq.gz > %o:$TMPDIR/test2.out 2> %o:$TMPDIR/test2.err"},
    {param=>"-g $INDIR/mxtest_2.fastq $INDIR/mxtest_1.fastq $INDIR/mxtest_3.fastq -o $TMPDIR/mxout_%_1.fq -o $TMPDIR/mxout_%_2.fq > %o:$TMPDIR/test3.out 2> %o:$TMPDIR/test3.err"},
    {param=>"-H -v ' ' -l $INDIR/master-barcodes.txt $INDIR/mxtest-h_1.fastq $INDIR/mxtest-h_2.fastq -o $TMPDIR/mxout_%_1.fq -o $TMPDIR/mxout_%_2.fq > %o:$TMPDIR/test4.out 2> %o:$TMPDIR/test4.err"},
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
