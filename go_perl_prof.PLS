#!/usr/bin/perl

# go_perl_prof.PLS
#
# Cared for by Albert Vilella <>
#
# Based on specification of the gene_association format is defined at:
#   http://www.geneontology.org/GO.annotation.html#file
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

go_perl_prof.PLS - DESCRIPTION 

=head1 SYNOPSIS

perl go_perl_prof.PLS  \
    -g /home/avb/wallace/go

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
use Getopt::Long;
use GO::Parser;
use File::Path;

my ($inputfile,
    $genelist,
    $outdir,
    $godir,
   );

$godir = "/home/avb/wallace/go";

GetOptions(
	   'i|input|inputfile:s' => \$inputfile,
	   'o|outdir:s' => \$outdir,
	   'g|godir:s' => \$godir,
	   'l|list|genelist:s' => \$genelist,
          );

my ($parser,$graph,$term,$ancestor_terms);
print "Loading ontology (might take a couple of minutes)...";
$parser = new GO::Parser({handler=>'obj'});
$parser->parse("$godir/gene_ontology.obo");
$graph = $parser->handler->graph;

1;;


