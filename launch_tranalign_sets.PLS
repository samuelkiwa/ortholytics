#!/usr/local/bin/perl
# launch_tranalign_sets.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

launch_tranalign_sets.PLS - DESCRIPTION 

=head1 SYNOPSIS

launch_tranalign_sets.PLS -dir /my/dir -aligndir probcons -outdir probcons_aa

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
use Bio::Factory::EMBOSS;

BEGIN {$ENV{EMBOSSDIR} = '/home/avb/bioperl-runnables/EMBOSS-2.9.0/emboss'; }

my ($factory, $tranalign, $result, $emboss_present,
    $nucl_file, $aa_aln, $tranalignoutfile, 
    $guessed_format, $aligndir) = "";
my ($dir, $outdir, $subdir, $file, $input, $seq, 
    $aln, @params, $factory, $out) = "";

GetOptions(
	   'd|dir:s' => \$dir,
	   'a|aligndir:s' => \$aligndir,
	   'o|outdir:s' => \$outdir,
          );

# get an EMBOSS application object from the EMBOSS factory
$factory = new Bio::Factory::EMBOSS;
$tranalign = $factory->program('tranalign');

# $emboss_present = $tranalign->executable();

# unless ($emboss_present) {
#     warn "emboss not found.\n";
#     exit(0);
# }

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
            while (( defined($aa_aln = readdir(ALIGNDIR)) )) {
                if ($aa_aln =~ /(\S+)\.afa$/) {
                        # running the application tranalign
                    $tranalignoutfile = "$dir/$subdir/$outdir/$nucl_file.nucl";
                    unless (-d "$dir/$subdir/$outdir") {
                        mkdir "$dir/$subdir/$outdir" or die "couldnt create directory: $!";
                    }
                    print "Processing\n$dir/$subdir/$nucl_file \n$dir/$subdir/$aligndir/$aa_aln ...\n";
                    $tranalign->run
                        (
                         {
                          '-asequence' => "$dir/$subdir/$nucl_file",
                          '-bsequence'    => "$dir/$subdir/$aligndir/$aa_aln",
                          '-outseq'    => $tranalignoutfile,
                         }
                        );
                    $aa_aln = "";
                    $nucl_file = "";
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


