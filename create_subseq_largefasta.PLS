#!/usr/local/bin/perl

# create_subseq_largefasta.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

create_subseq_largefasta.PLS - DESCRIPTION 

=head1 SYNOPSIS

perl create_subseq_largefasta.PLS -i sequence.fasta [--subseq 11231:34534]

=head1 DESCRIPTION

This script will read a largefasta file, prompt its id and length, and
if subseq option is defined, write down another fasta file with the
desired subseq.

=head1 AUTHOR - Albert Vilella

Email 

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut


# Let the code begin...

use strict;
use Bio::SeqIO;
use Bio::Seq::LargePrimarySeq;
use Getopt::Long;

my ($infile,$coordsfile,$subseq,$slice,$split,$quiet);

GetOptions(
	   'i|infile:s' => \$infile,
           'c|coords:s' => \$coordsfile,
	   's|sub|subseq:s' => \$subseq,
	   'slice:s' => \$slice,
	   'split' => \$split,
	   'quiet' => \$quiet,
          );

my ($start,$end) = split(/\:/, $subseq) if ($subseq);

my $in  = Bio::SeqIO->new
    (
     -file => "$infile" , 
     -format => 'largefasta'
    );

my %entries;
if ($coordsfile) {
    open COORDSFILE, "$coordsfile" or die "error opening $coordsfile: $!";
    while (<COORDSFILE>) {
        $_ =~ /(\d+),(\d+),(\d+),([\+\-]),(\w+)/;
        $entries{$1}{_start} = $2;
        $entries{$1}{_end} = $3;
        $entries{$1}{_strand} = $4;
        $entries{$1}{_contig} = $5;
    }
}

print "Reading sequences from file $infile (may take a while)...\n";
my $out = Bio::SeqIO->new
    (
     -file => ">$infile.fasta",
     -format => 'largefasta'
    ) unless ($split);
while ( my $seq = $in->next_seq() ) {
    my $id = $seq->id;
    print "Sequence ",$seq->id, " length ", $seq->length, "\n";
    if ($subseq) {
#         my $sub_seq = $seq->subseq($start,$end);
        my $out_seq = Bio::Seq::LargePrimarySeq->new
            (
             -id => $id,
             -alphabet => 'dna',
             -seq => $seq->subseq($start,$end),
            );
        print "Writing output file $infile.$id.$start.$end.fasta \n";
        my $out = Bio::SeqIO->new
            (
             -file => ">$infile.$id.$start.$end.fasta",
             -format => 'largefasta'
            );
        $out->write_seq($out_seq);
    } elsif (defined $slice) {
        my $length = $seq->length;
        my %slices;
        my $remains = $length;
        while ($remains > 0) {
            $remains = $remains - $slice;
            $remains = sprintf("%09d",$remains);
            $slices{$remains}{_start} = $remains - 1;
            $slices{$remains}{_end} = $slice + $remains;
            $slices{$remains}{_start} = 1 if ($slices{$remains}{_start} < 0);
        }
        foreach my $slice (sort keys %slices) {
            my $start = $slices{$slice}{_start};
            my $end = $slices{$slice}{_end};
            my $out_seq = Bio::Seq::LargePrimarySeq->new
                (
                 -id => $id,
                 -alphabet => 'dna',
                 -seq => $seq->subseq($start,$end),
                );
            $start = sprintf("%09d",$start);
            $end = sprintf("%09d",$end);
            print "Writing output file $infile.$start.$end.fasta \n";
            my $out = Bio::SeqIO->new
                (
                 -file => ">$infile.$start.$end.fasta",
                 -format => 'largefasta'
                );
            $out->write_seq($out_seq);
        }
    } elsif (defined %entries) {
        foreach my $entry (keys %entries) {
            if (defined $entries{$entry}{_contig}) {
                if ($id eq $entries{$entry}{_contig}) {
                    my $myseq;
                    $myseq = $seq->subseq($entries{$entry}{_start},$entries{$entry}{_end});
                    my $seqobj = Bio::PrimarySeq->new
                        ( -seq => $myseq,
                          -id  => $entry,
                          -moltype => 'dna',
                        );
                    $seqobj = $seqobj->revcom if ($entries{$entry}{_strand} ne "+");
                    $myseq = $seqobj->seq;
                    my $out_seq = Bio::Seq::LargePrimarySeq->new
                        (
                         -id => $entry,
                         -alphabet => 'dna',
                         -seq => $myseq,
                        );
                    print "$entry\n" unless($quiet);
                    $out->write_seq($out_seq);
                }
            }
        }
    } elsif ($split) {
        my $out = Bio::SeqIO->new
            (
             -file => ">$infile.$id.fasta",
             -format => 'largefasta'
            );
        $out->write_seq($seq);
        print "Out seq $infile.$id.fasta\n";
    }
}

1;;




