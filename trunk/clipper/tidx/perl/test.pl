use Text::Tidx;

my $cytoTab=Text::Tidx->new("/mnt/scratch/indexes/Annotations/cytoBand.txt");

print $cytoTab->lookup("chr1", 2300001), "\n";
print $cytoTab->query("chr1", 2300001), "\n";
