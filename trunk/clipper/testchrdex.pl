BEGIN {
	@INC = ('.', @INC);
}

use Chrdex;
use Data::Dumper;

#$x = Chrdex->new("chrdexannot.txt", chr=>2, beg=>5, end=>6, skip=>1);
#for ($i=0;$i<10000;++$i) {
#	$x->search(1, 153432255);
#}

my $ct=Chrdex->new("/mnt/scratch/indexes/Annotations/cytoBand.txt", force=>1);
print Dumper($ct->query(1, 8601118));

$x = Chrdex->new("chrdexannot.txt", chr=>2, beg=>5, end=>6, skip=>1);
print "---l1:\n";
print $x->search(1, 851184), "\n";
print "---l2:\n";
print $x->search(1, 869395), "\n";
print "---l2-3:\n";
print $x->search(1, 869395, 869936), "\n";
