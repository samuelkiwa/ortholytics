#!/usr/local/bin/perl
# launch_amap_sets.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

launch_amap_sets.PLS -

=head1 SYNOPSIS

nice perl launch_amap_sets.PLS -dir /my/dir -outdir amap -table 11 -email "myemail@address.com"

=head1 DESCRIPTION

launch_amap_sets will load the sequences in the subdirectories of
the speficied directory, translate them to aa, run amap, and give
flush the output alignment in the same subdirectory. 

Only works when the next-depth level subdir contains the files

-email "myemail@address.com" Sent an email when finished

=head1 AUTHOR - Albert Vilella

Email 

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut


# Let the code begin...

use strict;
use Getopt::Long;
use Bio::Tools::Run::Alignment::Amap;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Tools::GuessSeqFormat;
# use Bio::Tools::CodonTable;
# use Mail::Sendmail;

# BEGIN {$ENV{AMAPDIR} = '/home/avb/9_opl/amap/amap/'; }
BEGIN {$ENV{AMAPDIR} = '/home/avilella/9_opl/amap/amap/'; }


my ($dir, $outdir, $subdir, $file, $input, $aa_input, $seq, 
    @seq_array, $seq_array_ref, $aln, $pseq, $email,
    $codontable, $tableid, @params, $factory, $out, 
#     $matrix
   ) = "";

GetOptions(
 	   'table|tableid:s'=> \$tableid,
	   'dir:s' => \$dir,
	   'o|outdir:s' => \$outdir,
	   'e|email:s' => \$email,

          );

# $codontable = Bio::Tools::CodonTable->new( -id => $tableid);

opendir DIR, $dir or die "couldnt find $dir:$!";
while (defined($subdir = readdir(DIR))) {
    next if ($subdir =~ /\./);
    next if ($subdir =~ /\.\./);
    #FIXME: use dir more politely
    opendir SUBDIR, "$dir/$subdir" or die "couldnt find subdir $subdir:$!";
    while (defined($file = readdir(SUBDIR))) {
        if ($file =~ /(\S+)\.fasta$/) {
            print "Aligning $subdir\n";
            unless (-d "$dir/$subdir/$outdir") {
                mkdir "$dir/$subdir/$outdir" or die "couldnt create directory: $!";
            }
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
#                 $pseq = $seq->translate(undef, undef, undef, undef, undef, $tableid);
                $pseq = $seq->translate(undef, undef, undef, $tableid, 'fullCDS', undef, undef);
                $aa_input->write_seq($pseq);
                push (@seq_array, $pseq);
            }

            # Build a amap alignment factory
            $factory = new Bio::Tools::Run::Alignment::Amap
                (
                 'iterative-refinement'  => '1000',
                 'consistency'   => '5',
                 'emissions' => '',
                 'verbose' => '',
                 'train'   => "$dir/$subdir/$outdir/train.params",
                );
            $factory->outfile_name("$dir/$subdir/$outdir/train.params");

            my $amap_present = $factory->executable();
            unless ($amap_present) {
                warn "amap program not found.\n";
                exit(0);
            }

            $seq_array_ref = \@seq_array;
            $aln = $factory->align($seq_array_ref);
            $aln = '';
            $factory = '';

            $factory = new Bio::Tools::Run::Alignment::Amap
                (
                 'iterative-refinement'  => '1000',
                 'consistency'   => '5',
                 'verbose' => '',
                 'annot'   => "$dir/$subdir/$outdir/annotation.scores",
                 'paramfile'   => "$dir/$subdir/$outdir/train.params",
#                  'matrixfile'   => "$matrix",
                );
            $factory->outfile_name("$dir/$subdir/$outdir/temp");
            $aln = $factory->align($seq_array_ref);

#             $out = Bio::AlignIO->new(-file=>Bio::Root::IO->catfile(">","$dir","$subdir","$outdir","$file.phylip"),
#                                      -format => 'phylip');
            $out = Bio::AlignIO->new(-file=>Bio::Root::IO->catfile(">","$dir","$subdir","$outdir","$file.afa"),
                                     -format => 'fasta');
            $out->write_aln($aln);
            @seq_array = ();
            unlink("$dir/$subdir/$outdir/temp");
        }
    }
}

my $message = "Output ready at $dir/whatever/$outdir";

if ($email) {
    my %mail = ( To      => "$email",
                 From    => 'launch_amap_sets@localhost',
                 Subject => "Sendmail batch: $0 job finished.",
                 Message => "$message"
               );

    sendmail(%mail) or die $Mail::Sendmail::error;
    print "Mail sent. Log says:\n", $Mail::Sendmail::log;
    print "\n";
}

# eval {
# #something;
# };
# $@ ? ok 1 : ok 0;

1;
