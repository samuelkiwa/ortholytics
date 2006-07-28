#!/usr/local/bin/perl

# rst_seqs_paml.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

rst_seqs_paml.PLS - DESCRIPTION 

=head1 SYNOPSIS

perl rst_seqs_paml.PLS
-mlcfile /my/paml/results/dir/mlc.txt
-dir /my/paml/results/dir

=head1 DESCRIPTION

rst_seqs_paml.PLS will load the results of a PAML mlc file and
directory (where rst file is found), and create a fasta file for each
of the reconstructed seqs.

=head1 AUTHOR - Albert Vilella

Email 

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut


# Let the code begin...

use strict;
use Getopt::Long;
use Bio::Tools::Phylo::PAML;
use Bio::SeqIO;

my ($mlcfile, $dir, $merged);

GetOptions(
 	   'mlcfile|mlc|m:s'=> \$mlcfile,
	   'dir|d:s' => \$dir,
	   'merged' => \$merged,
          );

my $parser = new Bio::Tools::Phylo::PAML
  (
   -file => "$mlcfile",
   -dir => "$dir",
  );

# get the first/next result; a Bio::Tools::Phylo::PAML::Result object,
# which isa Bio::SeqAnalysisResultI object.
print "Parsing result...\n";
my $result = $parser->next_result();

# my @seqs         = $result->get_seqs;
# my %input_params = $result->get_input_parameters;
# my @basfreq      = $result->get_codon_pos_basefreq;
# my $MLmatrix     = $result->get_MLmatrix; # get MaxLikelihood Matrix
# my $NGmatrix     = $result->get_NGmatrix; # get Nei-Gojoburi Matrix

# if -dir contains an rst file get list of
# Bio::PrimarySeq ancestral state reconstructions of the sequences
print "Parsing rst seqs...\n";
my @rsts = $result->get_rst_seqs;
my $num = scalar @rsts;
print "$num rsts found...\n";

my $counter = 0;
my $rst_seq;
foreach my $rst (sort @rsts) {
    my $id = $rst->display_id;
    # Make it like evolver
    $id =~ s/\#//g;
    $id = $rst->display_id("$id");
    if (0 == $counter || !$merged) {
        $rst_seq = Bio::SeqIO->new
            (-file=>Bio::Root::IO->catfile(">","$dir","$id.fasta"),
             -format=>'fasta',
            );
    }
    my $desc = $rst->desc();
    $rst->desc('');
    print "Writing $id - $desc\n";
    $rst_seq->write_seq($rst);
    $counter++;
}

1;
