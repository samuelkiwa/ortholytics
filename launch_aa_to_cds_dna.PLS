#!/usr/local/bin/perl
# launch_aa_to_cds_dna_sets.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

launch_aa_to_cds_dna_sets.PLS - DESCRIPTION 

=head1 SYNOPSIS

launch_aa_to_cds_dna_sets.PLS -dir /my/dir -aligndir probcons -outdir probcons_aa [-removegaps] [-phylip]

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
use Bio::Root::IO;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Tools::GuessSeqFormat;
use Bio::Align::Utilities qw(aa_to_dna_aln);

my ($factory, $tranalign, $result, $emboss_present,
    $dna_aln, $dna_aln_nogaps, $aa_aln, $aa_aln_file, $outfile, $nucl_file, $input,
    $aa_input, $phylip,
    $guessed_format, $aligndir) = "";
my ($dir, $outdir, $subdir, $removegaps, $file, $seq, %seqs,
    $out) = "";

GetOptions(
	   'd|dir:s' => \$dir,
	   'a|aligndir:s' => \$aligndir,
	   'o|outdir:s' => \$outdir,
           'p|phy|phylip:s' => \$phylip,
           'r|rem|remove|removegaps' => \$removegaps,
          );

opendir DIR, $dir or die "couldnt find $dir:$!";
while (defined($subdir = readdir(DIR))) {
    next if ($subdir =~ /\./);
    next if ($subdir =~ /\.\./);
    #FIXME: use dir more politely
    opendir SUBDIR, "$dir/$subdir" or die "couldnt find subdir $subdir:$!";
    while (defined($nucl_file = readdir(SUBDIR))) {
        #FIXME: hardcoded fasta format
        if ($nucl_file =~ /(\S+)\.fasta$/) {
            opendir ALIGNDIR, "$dir/$subdir/$aligndir" or die "couldnt find subdir $aligndir:$!";
            while (( defined($aa_aln_file = readdir(ALIGNDIR)) )) {
                if ($aa_aln_file =~ /(\S+)\.afa$/) {
                    #Create directory
                    unless (-d "$dir/$subdir/$outdir") {
                        mkdir "$dir/$subdir/$outdir" or die "couldnt create directory: $!";
                    }

                    print "Processing\n$dir/$subdir/$nucl_file \n$dir/$subdir/$aligndir/$aa_aln_file ...\n";
                    #Load sequences
                    my $guessed_format = new Bio::Tools::GuessSeqFormat
                        (-file=>Bio::Root::IO->catfile("$dir","$subdir","$nucl_file")
                         #-verbose=> $verbose;
                        )->guess;
                    $input = Bio::SeqIO->new
                        (
                         -file=>Bio::Root::IO->catfile("$dir","$subdir","$nucl_file"),
                         -format=>$guessed_format,
                        );

                    my $aa_guessed_format = new Bio::Tools::GuessSeqFormat
                        (-file=>Bio::Root::IO->catfile("$dir","$subdir","$aligndir","$aa_aln_file")
                         #-verbose=> $verbose;
                        )->guess;
                    $aa_input = Bio::AlignIO->new
                        (
                         -file=>Bio::Root::IO->catfile("$dir","$subdir","$aligndir","$aa_aln_file"),
                         -format=>$aa_guessed_format,
                        );

                    $aa_aln = $aa_input->next_aln;

                    while($seq = $input->next_seq ) {
                        $seqs{$seq->id} = $seq;
                    }

                    $dna_aln = aa_to_dna_aln($aa_aln,\%seqs);
                    $dna_aln->set_displayname_flat(1);
                    if ($removegaps) {
                        $dna_aln_nogaps = $dna_aln->remove_gaps();
                    }

                    if ($phylip) {
                        $out = new Bio::AlignIO
                            (
                             -file=>Bio::Root::IO->catfile(">", "$dir","$subdir","$outdir","$aa_aln_file.phy"),
                             -format=> 'phylip'
                             -idlength => '50',
                             -interleaved => 0,
                            );
                    } else {
                        $out = new Bio::AlignIO
                            (
                             -file=>Bio::Root::IO->catfile(">", "$dir","$subdir","$outdir","$aa_aln_file.dna.fasta"),
                             -format=> 'fasta',
                            );
                    }

                    if ($removegaps) {
                        $out->write_aln($dna_aln_nogaps);
                    } else {
                        $out->write_aln($dna_aln);
                    }
                    %seqs = '';
                } 
            }
        }
#         $subdir = "";
    }
}

close DIR;
close SUBDIR;
close ALIGNDIR;


1;


