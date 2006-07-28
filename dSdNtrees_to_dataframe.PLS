#!/usr/local/bin/perl
# dSdNtrees_to_dataframe.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

dSdNtrees_to_dataframe.PLS - DESCRIPTION 

=head1 SYNOPSIS

perl \
dSdNtrees_to_dataframe.PLS \
-i \
/home/avb/9_opl/hyphy_latest/data/HYPHYMP0001/res.dSdNtrees.txt \
--tag \
teto

=head1 DESCRIPTION

Describe the object here

=head1 AUTHOR - Albert Vilella

Email 

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut


# Let the code begin...

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
my @entries;
foreach my $nodesyn (@all_nodes_syn) {
    # skip root node
    unless ($nodesyn->ancestor) { 
        $j++;
        next;
    }
    my $nodenon = $all_nodes_non[$j];
    my $asyn;
    $asyn = $nodesyn->ancestor->id;
    my $internalasyn = $nodesyn->ancestor->internal_id;
    unless ($asyn) {$asyn = $internalasyn};
    my $idsyn = $nodesyn->id;
    my $internalidsyn = $nodesyn->internal_id;
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
    } else {print "$tag;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;\n";}
    #   print "$tag;ancetorid..id;branch;t;S;N;dN/dS;dN;dS;S*dS;N*dN;\n";
    push @entries,"$tag;$id..$asyn;NaN;$asyn..$internalidsyn;NaN;NaN;NaN;$omega;$blnon;$blsyn;NaN;NaN;\n";
    $j++;
}

for my $entry (sort @entries) {
    print "$entry";
}

1;

