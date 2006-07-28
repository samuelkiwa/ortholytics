#!/usr/local/bin/perl

# gff2locus_GI_csv.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

gff2locus_GI_csv.PLS - DESCRIPTION 

=head1 SYNOPSIS

perl /home/avb/wallace/prokarya_archaea/scripts/gff2locus_GI_csv.PLS NC_006624.gff > NC_006624_locus_GI.csv

=head1 DESCRIPTION

Print a comma-separated list of locus_ids and NCBI-GI from a gff file

=head1 AUTHOR - Albert Vilella

Email 

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut


# Let the code begin...

use strict;

my $fh;
my $file = shift @ARGV;
open($fh, "$file") || die "$!";

my $lasttag = 0;
while(<$fh> ){
     if( /locus_tag\=(\w+)\;.+db_xref\=GI:(\d+)/ ) {
         if ($lasttag == 0) {
             print "$1,$2\n";
             $lasttag = 1;
         } else {
             $lasttag = 0;
         }
     }
}

close $fh;

1;
