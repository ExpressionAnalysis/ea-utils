#BEGIN {
#    @INC=("blib/lib",@INC);
#    $ENV{LD_LIBRARY_PATH}="./blib/arch/auto/Text/Tidx/SWIG";
#};

use Text::Tidx;

Text::Tidx::build("chr20.txt", "\t", 0, 3, 4, 0, '#');

$x = Text::Tidx->new("chr20.txt");

#$x=Text::Tidx::read("chr20.txt");

$v = $x->lookup("chr20", 251858, "^");

print "found: $v\n";
