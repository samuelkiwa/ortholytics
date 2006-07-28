#!/usr/local/bin/perl
use Getopt::Long;

my ($inputfile,$yesnoopt);

GetOptions(
	   'i|input|inputfile:s' => \$inputfile,
           'yesnoopt' => \$yesnoopt,
          );
use Bio::TreeIO;
my $in = new Bio::TreeIO(-file => $inputfile, -format => 'newick');
my $out = new Bio::TreeIO(-file => ">$inputfile.tabtree", -format => 'tabtree');

while( my $tree = $in->next_tree ) {
    $out->write_tree($tree);
}
1;;
