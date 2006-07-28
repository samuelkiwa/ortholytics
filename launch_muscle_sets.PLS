#!/usr/local/bin/perl
# launch_muscle_sets.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

launch_muscle_sets.PLS -  

=head1 SYNOPSIS

perl launch_muscle_sets.PLS 
-dir /my/dir 
-aligndir probcons 
-outdir muscle

=head1 DESCRIPTION

launch_muscle_sets will load the sequences in the subdirectories of
the speficied directory, translate them to aa, run muscle, and give
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

my ($dir, $outdir, $subdir, $file, $input, $aa_input, $seq, 
    @seq_array, $seq_array_ref, $aln, $pseq,
    $codontable, $tableid, @params, $factory, $out) = "";

GetOptions(
 	   'table|tableid:s'=> \$tableid,
	   'dir:s' => \$dir,
	   'o|outdir:s' => \$outdir,
          );

opendir DIR, $dir or die "couldnt find $dir:$!";
while (defined($subdir = readdir(DIR))) {
    next if ($subdir =~ /\./);
    next if ($subdir =~ /\.\./);
    #FIXME: use dir more politely
    opendir SUBDIR, "$dir/$subdir" or die "couldnt find subdir $subdir:$!";
    while (defined($file = readdir(SUBDIR))) {
        if ($file =~ /(\S+)\.fasta$/) {
            #Load sequences
            my $guessed_format = new Bio::Tools::GuessSeqFormat
                (-file=>Bio::Root::IO->catfile("$dir","$subdir","$file")
             #-verbose=> $verbose;
                )->guess;
            $input = Bio::SeqIO->new
                (-file=>Bio::Root::IO->catfile("$dir","$subdir","$file"),
                 -format=>$guessed_format,
                );
            $aa_input = Bio::SeqIO->new
                (-file=>Bio::Root::IO->catfile(">","$dir","$subdir","$file.aa"),
                 -format=>$guessed_format,
                );
            while($seq = $input->next_seq()){
                $pseq = $seq->translate(undef, undef, undef, $tableid, 'fullCDS', undef, undef);
                $aa_input->write_seq($pseq);
                push (@seq_array, $pseq);
                my $id = $pseq->display_id();
                my $length = $pseq->length();
                print ("$subdir;$id;$length\n");
            }

            $seq_array_ref = \@seq_array;
            $factory = '';
            unless (-d "$dir/$subdir/$outdir") {
                mkdir "$dir/$subdir/$outdir" or die "couldnt create directory: $!";
            }
            $factory = new Bio::Tools::Run::Alignment::Muscle
                (
                 'maxiters'  => '32000',
                 'tree1'   => "$dir/$subdir/$outdir/$subdir.tree.dnd",
#                  'scorefile'   => "$dir/$subdir/$outdir/$subdir.scores",
                );
            $aln = $factory->align($seq_array_ref);
            unless (-d "$dir/$subdir/$outdir") {
                mkdir "$dir/$subdir/$outdir" or die "couldnt create directory: $!";
            }
#             $out = Bio::AlignIO->new(-file=>Bio::Root::IO->catfile(">","$dir","$subdir","$outdir","$file.phylip"),
#                                      -format => 'phylip');
            $out = Bio::AlignIO->new(-file=>Bio::Root::IO->catfile(">","$dir","$subdir","$outdir","$file.afa"),
                                     -format => 'fasta');
            $out->write_aln($aln);
            @seq_array = ();
        }
    }
}

1;
