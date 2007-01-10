#!/usr/local/bin/perl

# Let the code begin...

use strict;
use Getopt::Long;
use Bio::Root::IO;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Tools::GuessSeqFormat;
use Bio::Align::Utilities qw(aa_to_dna_aln);
use Bio::Align::DNAStatistics;

my ($factory, $tranalign, $result, $emboss_present,
    $dna_aln, $dna_aln_nogaps, $aa_aln, $aa_aln_file, $outfile, $id_file, $input,
    $aa_input, $phylip, $comparedirs, $guessed_format, $debug) = "";
my ($dir, $outdir, $subdir, $removegaps, $file, $seq, %seqs,
    $out) = "";

GetOptions(
	   'd|dir:s' => \$dir,
	   'c|comparedirs:s' => \$comparedirs,
	   'o|outdir:s' => \$outdir,
           'p|phy|phylip:s' => \$phylip,
           'r|rem|remove|removegaps' => \$removegaps,
           'debug:s' => \$debug,
          );

my $stats;
foreach my $comparedir (split ":", $comparedirs) {
    $stats->{_comparedirs}{$comparedir} = 1;
}

opendir DIR, $dir or die "couldnt find $dir:$!";
my @subdirs = readdir(DIR);
my $count=0;
foreach my $subdir (sort @subdirs) {
    next if ($subdir =~ /\./);
    next if ($subdir =~ /\.\./);
    #FIXME: use dir more politely
    $count++; $count = sprintf("%04d", $count);
    next if ($count < $debug);
    print STDERR "$count -- $dir/$subdir\n";
    opendir SUBDIR, "$dir/$subdir" or die "couldnt find subdir $dir/$subdir:$!";
    my @aligndirs = readdir(SUBDIR);
    foreach my $aligndir (@aligndirs) {
        next unless (defined($stats->{_comparedirs}{$aligndir}));
        my $missing = 0;
        unless(opendir ALIGNDIR, "$dir/$subdir/$aligndir") {
            warn "\nmissing alignment $dir/$subdir/$aligndir\n";
            $missing = 1;
        }
        next if ($missing);
        while (( defined($aa_aln_file = readdir(ALIGNDIR)) )) {
            if ($aa_aln_file =~ /(\S+)\.afa$/) {
                #Create directory
                unless (-d "$dir/$subdir/$outdir") {
                    mkdir "$dir/$subdir/$outdir" or die "couldnt create directory: $!";
                }
                print STDERR "     -- $dir/$subdir/$id_file \n$dir/$subdir/$aligndir/$aa_aln_file ...\n";
                #Load sequences
                my $guessed_format = new Bio::Tools::GuessSeqFormat
                    (-file=>Bio::Root::IO->catfile("$dir","$subdir","$aligndir","$aa_aln_file")
                     #-verbose=> $verbose;
                    )->guess;
                $input = Bio::SeqIO->new
                    (
                     -file=>Bio::Root::IO->catfile("$dir","$subdir","$aligndir","$aa_aln_file"),
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
                unless ($aa_aln->is_flush) {
                    print STDERR "check $subdir/$aligndir\n";
                    next;
                }
                $stats->{_alndirs}{$subdir}{$aligndir}{_overall_percid} = $aa_aln->overall_percentage_identity;
                $stats->{_alndirs}{$subdir}{$aligndir}{_average_percid} = $aa_aln->average_percentage_identity;
                $stats->{_alndirs}{$subdir}{$aligndir}{_alnlength} = $aa_aln->length;
                $stats->{_alndirs}{$subdir}{$aligndir}{_patterns} = $aa_aln->remove_gaps->length;
                print $subdir, 
                ", ",
                $aligndir, 
                ", ", 
                $stats->{_alndirs}{$subdir}{$aligndir}{_overall_percid}, 
                ", ", 
                $stats->{_alndirs}{$subdir}{$aligndir}{_average_percid}, 
                ", ", 
                $stats->{_alndirs}{$subdir}{$aligndir}{_patterns}, 
                ", ", 
                $stats->{_alndirs}{$subdir}{$aligndir}{_alnlength}, "\n";
            }
        }
        1;
    }
    1;
}

close DIR;
close SUBDIR;
close ALIGNDIR;

# foreach my $subdir (keys %{$stats->{_alndirs}}) {
#     foreach my $aligndir (keys %{$stats->{_alndirs}{$subdir}}) {
#         print $subdir, 
#             ", ",
#             $aligndir, 
#             ", ", 
#             $stats->{_alndirs}{$subdir}{$aligndir}{_overall_percid}, 
#             ", ", 
#             $stats->{_alndirs}{$subdir}{$aligndir}{_average_percid}, "\n";
#             ", ", 
#             $stats->{_alndirs}{$subdir}{$aligndir}{_patterns}, "\n";
#             ", ", 
#             $stats->{_alndirs}{$subdir}{$aligndir}{_alnlength}, "\n";
#     }
# }

1;


