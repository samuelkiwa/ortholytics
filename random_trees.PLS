#!/usr/local/bin/perl
# random_trees.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

random_trees.PLS - DESCRIPTION 

=head1 SYNOPSIS

perl random_trees.PLS \
 -treestring "((((Dmel:LENGTH,Dsim:LENGTH),Dyak:LENGTH),Dana:LENGTH),Dpse:LENGTH)"

=head1 DESCRIPTION

Describe the object here

	   'string|treestring:s' => \$treestring,
	   'm|maxcount|i|iter|iteration|iterations:s' => \$maxcount,
           'r|range|:s'=> \$range,
           'o|out|output|outfile:s'=> \$outfile,

=head1 AUTHOR - Albert Vilella

Email 

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut


# Let the code begin...

use strict;
use Getopt::Long;

my ($taxa,$outfilename,$maxcount,$treestring,$tree,$range);

$tree = "output";
$maxcount = 10;
$range = 10;

GetOptions(
	   'string|treestring:s' => \$treestring,
	   'm|maxcount|i|iter|iteration|iterations:s' => \$maxcount,
           'r|range|:s'=> \$range,
           'o|out|output|outfile:s'=> \$outfilename,
          );

my $outfile = $outfilename || "$tree.$maxcount.tree";

print "Generating $maxcount trees\n";
print "Range of distances (0" . ".."."$range)\n";
print "Creating $outfile\n";

open OUTFILE, ">$outfile" or die "error opening $outfile: $!";
while ($maxcount) {
    my $randomtree = $treestring;
    while($randomtree =~ s{LENGTH}{sprintf("%6f",rand(10));}gex) {
        print OUTFILE "$randomtree\n";
    }
    $maxcount--;
}

close OUTFILE;
1;
