#!/usr/local/bin/perl
use Bio::AlignIO;

$inputfilename = shift;
$in  = Bio::AlignIO->new(-file   => $inputfilename ,
                         -format => 'phylip');

my $counter=1;
while ( my $aln = $in->next_aln() ) {
    1;;
    $out = Bio::AlignIO->new
        (
         -file   => ">$inputfilename.$counter.fasta",
         -format => 'fasta',
        );
    $out->write_aln($aln);
    $counter++;
}
