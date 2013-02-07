use Tree::Suffix;

$tree = Tree::Suffix->new;

$tree->insert("apple\t32", "banana\t45", "appley\t38");

@pos = $tree->find("appl");

for (@pos) {
	print "pos: $_\n"
}
