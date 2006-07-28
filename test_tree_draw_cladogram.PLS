#!/usr/local/bin/perl
use strict;
use Bio::Tree::Draw::Cladogram;
use Bio::TreeIO;
use Getopt::Long;

my ($input);

GetOptions(
	   'i|input:s' => \$input,
          );

my $treeio = new Bio::TreeIO('-format' => 'newick',
#                              '-file'   => '/home/avb/wallace/prokarya_archaea/pyrococcales/tmcales_concat/kadesh/data/hyphy01.trees');
                              '-file'   => $input);
#                              '-file'   => '/home/avb/wallace/prokarya_archaea/methanos/mbia_concat/8.7.HKY85.2.dnd');
my $t1 = $treeio->next_tree;
my $t2 = $treeio->next_tree;

my $tout = new Bio::TreeIO('-format' => 'newick',
                           -fh   => \*STDOUT,
                          );

# my $obj1 = new Bio::Tree::Draw::Cladogram(-tree => $t1);
# $obj1->print(-file => "$input.cladogram1.eps");

if ($t2) {
#     my $obj1_2 = new Bio::Tree::Draw::Cladogram(-tree => $t2);
#     $obj1->print(-file => "$input.cladogram2.eps");
    my $obj2 = new Bio::Tree::Draw::Cladogram
        (-tree => $t1,
         -second => $t2,
         -compact => 1,
         -size => 11,
#         -column => 10,
#         -xstep => 10,
#         -colors => 1,
        );
    $obj2->print(-file => "$input.tanglegram.eps");
}

1;
