#!/usr/local/bin/perl
use strict;
use Bio::TreeIO;
use Getopt::Long;

my ($inputfile,$tag);

GetOptions(
	   'i|input|inputfile:s' => \$inputfile,
	   't|tag:s' => \$tag,
          );

my $treeio = new Bio::TreeIO
    (
     '-format' => 'nhx',
     '-file'   => $inputfile,
    );

my $treesyn = $treeio->next_tree;
my $treenon = $treeio->next_tree;

my $total_length_syn = $treesyn->total_branch_length;
my $total_length_non = $treenon->total_branch_length;
my @all_nodes_syn = $treesyn->get_nodes;
my @all_nodes_non = $treenon->get_nodes;
my $num = sprintf("%02d",1);
my $j=0;
foreach my $nodesyn (@all_nodes_syn) {
    # skip root node
    unless ($nodesyn->ancestor) { 
        $j++;
        next;
    }
    my $nodenon = $all_nodes_non[$j];
    my $asyn;
    $asyn = $nodesyn->ancestor->id;
    unless ($asyn) {$asyn = $nodesyn->ancestor->internal_id;}
    my $idsyn = $nodesyn->id;
    my $idnon = $nodenon->id;
    my $blsyn = $nodesyn->branch_length;
    my $blnon = $nodenon->branch_length;
    my $omega = $blnon/$blsyn;
    $blsyn = sprintf("%.4f",$blsyn);
    $blnon = sprintf("%.4f",$blnon);
    $omega = sprintf("%.4f",$omega);
    my $id;
    if ($idsyn) {
        $id = $idsyn;
    } elsif (defined $blsyn) {
        $id = "Node$num";
        $num++;
    } else {print "error ; NaN\n";}
    #   print "$tag;ancetorid..id;branch;t;S;N;dN/dS;dN;dS;S*dS;N*dN;\n";
    print "$tag;$asyn..$id;nan;nan;nan;nan;$omega;$blnon;$blsyn;nan;nan;\n";
}

1;;

# my @leaf_nodes = $tree->get_leaf_nodes;
# my $num_internal = scalar(@all_nodes) - scalar(@leaf_nodes);
