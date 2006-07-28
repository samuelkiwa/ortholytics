#!/usr/local/bin/perl
# rename_aln_seqs.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

rename_aln_seqs.PLS - DESCRIPTION 

=head1 SYNOPSIS

perl \
rename_aln_seqs.PLS \
-i \
seq.fasta \
-namelist \
Dmel:Dsim:Dsec:Dyak:Dere:Dana:Dpse:Dvir:Dmoj:Dgri

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
use Bio::Tools::GuessSeqFormat;
use Bio::SeqIO;
use Bio::Seq::LargePrimarySeq;
use Bio::AlignIO;

my ($inputfile,$namelist,$quiet);

GetOptions(
	   'i|input|inputfile:s' => \$inputfile,
	   'namelist:s' => \$namelist,
           'quiet' => \$quiet,
          );

my @names = split /\:/,$namelist if ($namelist);

my $guessed_format = new Bio::Tools::GuessSeqFormat
    (-file=> "$inputfile")->guess;
my $in = Bio::SeqIO->new
    (-file=>"$inputfile",
     -format=>$guessed_format,
    );

my %seqs;
while ( my $seq = $in->next_seq() ) {
    print STDERR "Sequence ",$seq->id, " length ", $seq->length, "\n" unless($quiet);
    my $id = $seq->id;
    $seqs{$id}{_seq} = $seq;
}

my $out = Bio::SeqIO->new
    (-file=>">$inputfile.renamed.$guessed_format",
     -format=>$guessed_format,
    );

my $count = 0;
foreach my $seqid (sort keys %seqs) {
    my $newid = "$names[$count]";
    print STDERR "$seqid to $newid\n" unless($quiet);
    my $seq = $seqs{$seqid}{_seq};
    my $len = $seq->length;
    my $out_seq = Bio::Seq::LargePrimarySeq->new
        (
         -id => "$newid",
         -alphabet => 'dna',
         -seq => $seq->subseq(1,$len),
        );
    $out->write_seq($out_seq);
    $count++;
}
print STDERR "Output $inputfile.renamed.$guessed_format\n" unless($quiet);
1;
