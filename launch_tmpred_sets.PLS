#!/usr/local/bin/perl
# launch_tmpred_sets.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

launch_tmpred_sets.PLS -  

=head1 SYNOPSIS

perl launch_tmpred_sets.PLS 
-dir /my/dir 

=head1 DESCRIPTION

launch_tmpred_sets will load the sequences in the subdirectories of
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
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Tools::GuessSeqFormat;
# use Bio::Tools::CodonTable;

my $tmpred = '/home/avb/9_opl/tmpred/tmpred/linux_intel/tmpred';
my $tmpreddir = '/home/avb/9_opl/tmpred/tmpred/linux_intel/';

my ($dir, $outdir, $subdir, $file, $aa_input, $seq, 
    @seq_array, $seq_array_ref, $aln, $pseq,
    $codontable, $tableid, @params, $factory, $aa_file, $out) = "";

GetOptions(
	   'dir:s' => \$dir,
          );

opendir DIR, $dir or die "couldnt find $dir:$!";
while (defined($subdir = readdir(DIR))) {
    next if ($subdir =~ /\./);
    next if ($subdir =~ /\.\./);
    #FIXME: use dir more politely
    opendir SUBDIR, "$dir/$subdir" or die "couldnt find subdir $subdir:$!";
    while (defined($file = readdir(SUBDIR))) {
        if ($file =~ /(\S+)\.fasta.aa$/) {
            #Load sequences
            my $guessed_format = new Bio::Tools::GuessSeqFormat
                (-file=>Bio::Root::IO->catfile("$dir","$subdir","$file")
             #-verbose=> $verbose;
                )->guess;
            $aa_input = Bio::SeqIO->new
                (-file=>Bio::Root::IO->catfile("$dir","$subdir","$file"),
                 -format=>$guessed_format,
                );
            my $seq;
            while ($seq = $aa_input->next_seq) {
                my $id = $seq->display_id;
                $aa_file = Bio::SeqIO->new
                    (-file=>Bio::Root::IO->catfile(">$dir","$subdir","$id.fasta.aa"),
                     -format=>$guessed_format,
                    );
                $aa_file->write_seq($seq);
                my $filename = "$dir/$subdir/$id.fasta.aa";
#                 my $result = `$tmpred < "$filename\n\n17\n33\n\n\n" > "$filename.out"`;
                chdir $tmpreddir; 
                system("echo -e \"$filename\n\n17\n33\n\n\n\" | ./tmpred $filename > $filename.out");
#                 my $result = system("$tmpred $filename\n\n17\n33\n\n\n > $filename.out");
                1;
            }
        }
    }
}

1;
