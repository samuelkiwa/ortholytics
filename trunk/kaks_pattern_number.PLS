#!/usr/local/bin/perl
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
