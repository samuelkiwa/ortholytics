#!/usr/local/bin/perl
# launch_tcoffee_sets.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

launch_tcoffee_sets.PLS -  

=head1 SYNOPSIS

perl launch_tcoffee_sets.PLS -dir /my/dir

=head1 DESCRIPTION

launch_tcoffee_sets will load the sequences in the subdirectories of
the speficied directory, translate them to aa, run tcoffee, and give
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
use Bio::Tools::Run::Alignment::Tcoffee;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Tools::GuessSeqFormat;
# use Bio::Tools::CodonTable;

# BEGIN {$ENV{TCOFFEEDIR} = '/home/avb/9_opl/tcoffee/tcoffee3.52_src/'; }
BEGIN {$ENV{TCOFFEEDIR} = '/home/avb/9_opl/T-COFFEE_distribution_Version_1.37/bin/'; }


my ($dir, $outdir, $subdir, $file, $input, $seq, 
    @seq_array, $seq_array_ref, $aln, $aa_seq,
    $codontable, $tableid, @params, $factory, $out) = "";

GetOptions(
 	   'table|tableid:s'=> \$tableid,
	   'dir:s' => \$dir,
	   'o|outdir:s' => \$outdir,
          );

# $codontable = Bio::Tools::CodonTable->new( -id => $tableid);

# Build a tcoffee alignment factory
@params = ();
$factory = new Bio::Tools::Run::Alignment::Tcoffee(@params);

my $tcoffee_present = $factory->executable();
unless ($tcoffee_present) {
    warn "tcoffee program not found.\n";
    exit(0);
}

opendir DIR, $dir or die "couldnt find $dir:$!";
while (defined($subdir = readdir(DIR))) {
    next if ($subdir =~ /\./);
    next if ($subdir =~ /\.\./);
    #FIXME: use dir more politely
    opendir SUBDIR, "$dir/$subdir" or die "couldnt find subdir $subdir:$!";
    while (defined($file = readdir(SUBDIR))) {
        if ($file =~ /(\S+)\.(\S+)$/) {

            #Load sequences
            my $guessed_format = new Bio::Tools::GuessSeqFormat
                (-file=>Bio::Root::IO->catfile("$dir","$subdir","$file")
             #-verbose=> $verbose;
                )->guess;
            $input = Bio::SeqIO->new
                (-file=>Bio::Root::IO->catfile("$dir","$subdir","$file"),
                 -format=>$guessed_format,
                );
            while($seq = $input->next_seq()){
                $aa_seq = $seq->translate(undef, undef, undef, undef, $tableid);
                push (@seq_array, $aa_seq);
            }
            $seq_array_ref = \@seq_array;
            $aln = $factory->align($seq_array_ref);
            unless (-d "$dir/$subdir/$outdir") {
                mkdir "$dir/$subdir/$outdir" or die "couldnt create directory: $!";
            }
            $out = Bio::AlignIO->new(-file=>Bio::Root::IO->catfile(">","$dir","$subdir","$outdir","$file.afa"),
                                     -format => 'fasta');
            $out->write_aln($aln);
            @seq_array = ();
        }
    }
}

# eval {
# #something;
# };
# $@ ? ok 1 : ok 0;

1;
