use Getopt::Long;
GetOptions("count=i");
my ($f1, $f2) = @ARGV;
open IN1, $f1;
open IN2, $f2;
while($i1=<IN1>) {
    $s1=<IN1>;
    $c1=<IN1>;
    $q1=<IN1>;
    $i2=<IN2>;
    $s2=<IN2>;
    $c2=<IN2>;
    $q2=<IN2>;
    print $i1, $s1, $c1, $q1, $i2, $s2, $c2, $q2;
    $c+=2;
    last if $c >= $opt_count;
}

