#!/usr/local/bin/perl
# kaks_pattern_number.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

kaks_pattern_number.PLS - DESCRIPTION 

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

This script will take a CDS MSA and count the number of codons without
gaps in any sequence, compared to the total number of codons to the
alignment

=head1 AUTHOR - Albert Vilella

Email 

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut


# Let the code begin...
use Bio::AlignIO;
use Bio::Align::DNAStatistics;
use strict;
use Getopt::Long;
use File::Basename;

my ($inputfile,$yesnoopt);
GetOptions(
	   'i|input|inputfile:s' => \$inputfile,
           'yesnoopt' => \$yesnoopt,
          );
my ($infilebase,$path,$type) = fileparse($inputfile,qr{\..+});
my $stats = new Bio::Align::DNAStatistics;
my $in = new Bio::AlignIO(-format => 'phylip',
                          -file   => $inputfile);
my $alnobj = $in->next_aln;
my $total_codons = ($alnobj->length)/3;
my $patterns = $stats->kaks_pattern_number($alnobj);
print "$infilebase, $patterns, $total_codons\n";
1;
