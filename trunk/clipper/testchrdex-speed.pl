use Benchmark qw(:all);

BEGIN {
        @INC = ('.', @INC);
}

use Chrdex;

$x = Chrdex->new("chrdexannot.txt", chr=>2, beg=>5, end=>6, skip=>1);
$r = Chrdex->new("chrdexannot.txt", chr=>2, beg=>5, end=>6, skip=>1, index_path=>"chrdexbyref", byref=>1);

cmpthese( -1, { 
	ram =>sub{$x->query(1, 869395)}, 
	ram_range =>sub{$x->query(1, 869395, 869936)}, 
	ref =>sub{$r->query(1, 869395)}, 
	ref_range =>sub{$r->query(1, 869395, 869936)},
	} 
);
