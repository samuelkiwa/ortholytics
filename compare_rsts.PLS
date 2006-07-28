#!/usr/local/bin/perl
# compare_rsts.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

compare_rsts.PLS - DESCRIPTION 

=head1 SYNOPSIS

perl compare_rsts.PLS \
-dir /path/to/seqfiles \

=head1 DESCRIPTION

This script will load all the fasta files in the subdirectories of a
given directory and compare the sequences in them against the
ancestral.seq file using the Kimura 2p method

=head1 AUTHOR - Albert Vilella

Email 

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut


# Let the code begin...

use strict;
use Getopt::Long;
use Bio::Root::IO;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Align::DNAStatistics;
use Bio::Matrix::Scoring;

my ($dir,$yesnoopt);

GetOptions(
	   'd|dir|inputdir:s' => \$dir,
           'yesnoopt' => \$yesnoopt,
          );

my ($subdir,$fasta_file,$fasta_input,$fasta);

my $guessed_format = new Bio::Tools::GuessSeqFormat
    (-file=>Bio::Root::IO->catfile("$dir","ancestral.seq")
    )->guess;
my $ancestral_input = Bio::AlignIO->new
    (
     -file=>Bio::Root::IO->catfile("$dir","ancestral.seq"),
     -format=>$guessed_format,
    );
my $ancestral = $ancestral_input->next_aln();
my @ancestral_seqs = $ancestral->each_seq();
my %ancestral_seqs;
my %node_alns;
foreach my $seq (@ancestral_seqs) {
    my $node_aln = new Bio::SimpleAlign();
    my $id = $seq->display_id;
    $ancestral_seqs{$id} = $seq->seq();
    my $anc_locseq = new Bio::LocatableSeq
        (-seq => $ancestral_seqs{$id},
         -id  => "true_ances_$id",
        ) if (defined($ancestral_seqs{$id}));
    $node_aln->add_seq($anc_locseq);
    $node_alns{$id} = $node_aln;
}
untie(@ancestral_seqs);


opendir DIR, $dir or die "couldnt find $dir:$!";
while (defined($subdir = readdir(DIR))) {
    next if ($subdir =~ /\./);
    next if ($subdir =~ /\.\./);
    next unless (-d "$dir/$subdir");
    opendir SUBDIR, "$dir/$subdir" or die "couldnt find $dir/$subdir:$!";
    while (( defined($fasta_file = readdir(SUBDIR)) )) {
        #Look for the alignments in the dir
        if ($fasta_file =~ /(\w+)(\-\w+)?\.fasta$/) {
            my $guessed_format = new Bio::Tools::GuessSeqFormat
                (-file=>Bio::Root::IO->catfile("$dir","$subdir","$fasta_file")
                )->guess;
            $fasta_input = Bio::SeqIO->new
                (
                 -file=>Bio::Root::IO->catfile("$dir","$subdir","$fasta_file"),
                 -format=>$guessed_format,
                );
            while($fasta = $fasta_input->next_seq()) {
                my $fasta_seq = $fasta->seq();
                my $id = $fasta->display_id;
                my $locseq = new Bio::LocatableSeq
                    (-seq => $fasta_seq,
                     -id  => "$subdir"."_"."$id",
                    );
                $id =~ s/\#//g;
                $node_alns{$id}->add_seq($locseq);
            }
        }
    }
}

1;;

foreach my $id (keys %node_alns) {
    next unless (defined $id);
    my $stats = new Bio::Align::DNAStatistics;
    my $aln = $node_alns{$id};
    my $out_aln = Bio::AlignIO->new
    (
     -file=>Bio::Root::IO->catfile(">$dir","$id.all.fasta"),
     -format=>'fasta',
    );
    $out_aln->write_aln($aln);
    $out_aln->close();
    my $k2p_dist = $stats->D_Kimura($aln);
    #     my ($k2p_dist,$k2p_var) = $stats->D_Kimura_variance($aln);
    #     print "k2p_var\n";
    #     print $k2p_var->print_matrix;
    print "k2p_dist\n";
    print $k2p_dist->print_matrix;
    my @distances = split /\s+/, $k2p_dist->print_matrix;
    my $avg_k2p = 0;
    my $numc = 0;
    foreach my $entry (@distances) {
        if ($entry =~ /[0-9]*\.[0-9]*/) {
            unless (0.0 >= $entry) {
                $avg_k2p += $entry;
                $numc++;
            }
        }
    }
    $avg_k2p = $avg_k2p/($numc);
    print "Avg_k2p = $avg_k2p\n\n";
}

1;
