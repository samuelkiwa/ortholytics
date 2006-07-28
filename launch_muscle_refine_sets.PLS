#!/usr/local/bin/perl
# launch_muscle_refine_sets.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

launch_muscle_refine_sets.PLS -  

=head1 SYNOPSIS

perl launch_muscle_refine_sets.PLS 
-dir /my/dir 
-aligndir probcons 
-outdir muscle_refine

=head1 DESCRIPTION

launch_muscle_refine_sets will load the sequences in the subdirectories of
the speficied directory, translate them to aa, run muscle_refine, and give
flush the output alignment in the same subdirectory. 

Only works when the next-depth level subdir contains the files

=head1 AUTHOR - Albert Vilella

Email 

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut

# Let the code begin...

use strict;
use Getopt::Long;
use Bio::Tools::Run::Alignment::Muscle;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Tools::GuessSeqFormat;
# use Bio::Tools::CodonTable;

# BEGIN {$ENV{MUSCLEDIR} = '/home/avb/9_opl/muscle/muscle3.52_src/'; }
BEGIN {$ENV{MUSCLEDIR} = '/home/avb/9_opl/muscle/muscle3.6_src/'; }
# BEGIN {$ENV{MUSCLEDIR} = '/home/avb/9_opl/muscle/muscle6_bin/'; }

my ($factory, $aa_aln, $aa_aln_file, $outfile, $input,
    $aa_input, $guessed_format, $aligndir) = "";
my ($dir, $outdir, $subdir, $file, $seq, %seqs,
    $out) = "";

GetOptions(
	   'd|dir:s' => \$dir,
	   'a|aligndir:s' => \$aligndir,
	   'o|outdir:s' => \$outdir,
          );

opendir DIR, $dir or die "couldnt find $dir:$!";
while (defined($subdir = readdir(DIR))) {
    next if ($subdir =~ /\./);
    next if ($subdir =~ /\.\./);
    opendir ALIGNDIR, "$dir/$subdir/$aligndir" or die "couldnt find subdir $dir/$subdir/$aligndir:$!";
    while (( defined($aa_aln_file = readdir(ALIGNDIR)) )) {
        if ($aa_aln_file =~ /(\S+)\.afa(\S+)$/) {
            #Create directory
            unless (-d "$dir/$subdir/$outdir") {
                mkdir "$dir/$subdir/$outdir" or die "couldnt create directory: $!";
            }

            print "Processing $dir/$subdir/$aligndir/$aa_aln_file ...\n";
            $aa_input = Bio::Root::IO->catfile
                (
                 "$dir","$subdir","$aligndir","$aa_aln_file"
                );

#             my $aa_guessed_format = new Bio::Tools::GuessSeqFormat
#                 (-file=>Bio::Root::IO->catfile("$dir","$subdir","$aligndir","$aa_aln_file")
#                  #-verbose=> $verbose;
#                 )->guess;
#             $aa_input = Bio::AlignIO->new
#                 (
#                  -file=>Bio::Root::IO->catfile("$dir","$subdir","$aligndir","$aa_aln_file"),
#                  -format=>$aa_guessed_format,
#                 );
#             $aa_aln = $aa_input->next_aln;

            # Build a muscle alignment factory
            $out = Bio::Root::IO->catfile
                (
                 "$dir","$subdir","$outdir","refined_$aa_aln_file"
                );
            $factory = new Bio::Tools::Run::Alignment::Muscle
                (
                 'maxiters'  => '32000',
                 'refine' => '',
#                  'in' => $aa_input,
                );
            $factory->outfile_name("$out");

            my $muscle_present = $factory->executable();
            unless ($muscle_present) {
                warn "muscle program not found.\n";
                exit(0);
            }
            $out = $factory->align($aa_input);
        }
    }
}

close DIR;
close SUBDIR;
close ALIGNDIR;

1;
