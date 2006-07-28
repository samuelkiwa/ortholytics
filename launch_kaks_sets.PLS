#!/usr/local/bin/perl
# launch_kaks_sets.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

launch_kaks_sets.PLS - DESCRIPTION 

=head1 SYNOPSIS

launch_kaks_sets.PLS -dir /my/dir -aligndir probcons -outdir probcons_aa

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
use Bio::Tools::Run::Phylo::PAML::Codeml;
use Bio::Tools::Run::Phylo::PAML::Yn00;

BEGIN {$ENV{PAMLDIR} = '/home/avb/9_opl/paml/paml3.14/src'; }

my ($factory, $kaks_factory, $result,
    $dna_aln_file, $dna_aln, $dna_input,
    $aligndir) = "";
my ($dir, $outdir, $subdir, $verbose, $file, $input, $seq, 
    $aln, @params, $factory, $out, $logfile) = "";
my  $kaks_prog = 'codeml';

GetOptions(
	   'd|dir:s' => \$dir,
	   'a|aligndir:s' => \$aligndir,
	   'o|outdir:s' => \$outdir,
           'v|verbose' => \$verbose,
	   'kaks:s'         => \$kaks_prog,
          );

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $logfile_time = 
    sprintf ("%04d%02d%02d_%02d%02d%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec);

$logfile = Bio::Root::IO->new
    (
     -file=>Bio::Root::IO->catfile
     (
      ">","$dir","logfile.$kaks_prog.$logfile_time.txt"
     ),
    );

printf ("Writing log to %s\n", $logfile->file);

$logfile->_print
    (
     "Started $logfile_time\n"
    );

$logfile->_print
    (
     "With options:\n -dir $dir\n -aligndir $aligndir \n -outdir $outdir \n -kaks $kaks_prog\n"
    );

opendir DIR, $dir or die "couldnt find $dir:$!";
while (defined($subdir = readdir(DIR))) {
    next if ($subdir =~ /\./);
    next if ($subdir =~ /\.\./);
    opendir ALIGNDIR, "$dir/$subdir/$aligndir" or die "couldnt find subdir $aligndir:$!";
    while (( defined($dna_aln_file = readdir(ALIGNDIR)) )) {
        if ($dna_aln_file =~ /(\S+)\.dna.fasta$/) {

            # running the kaks application
            $logfile->_print("Processing\n$dir/$subdir/$aligndir/$dna_aln_file ...\n");
            unless (-d "$dir/$subdir/$outdir") {
                mkdir "$dir/$subdir/$outdir" or die "couldnt create directory: $!";
            }

            my $aa_guessed_format = new Bio::Tools::GuessSeqFormat
                (-file=>Bio::Root::IO->catfile("$dir","$subdir","$aligndir","$dna_aln_file")
                )->guess;
            $dna_input = Bio::AlignIO->new
                (
                 -file=>Bio::Root::IO->catfile("$dir","$subdir","$aligndir","$dna_aln_file"),
                 -format=>$aa_guessed_format,
                );

            $dna_input->idlength(50);
            $dna_input->interleaved(0);
            $dna_aln = $dna_input->next_aln;
            my @each = $dna_aln->each_seq();

            if ( $kaks_prog =~ /yn00/i ) {
                $kaks_factory = Bio::Tools::Run::Phylo::PAML::Yn00->new(-verbose => $verbose);
            } elsif ( $kaks_prog =~ /codeml/i ) {
                # change the parameters here if you want to tweak your Codeml running!
                $kaks_factory = Bio::Tools::Run::Phylo::PAML::Codeml->new
                    (-verbose => $verbose,
                     -params => { 'runmode' => -2,
                                  'seqtype' => 1,
                                }
                    );
            }

            unless ( $kaks_factory->executable ) {
                warn("Could not find the executable for $kaks_prog, define PAMLDIR");
                exit(0);
            }

            my $out_file = Bio::Root::IO->catfile
                (
                 "$dir","$subdir","$outdir","$kaks_prog.results.txt"
                );
            unless ((-f $out_file) && (-b $out_file)) {
                $out = Bio::Root::IO->new
                    (
                     -file=>Bio::Root::IO->catfile(">","$dir","$subdir","$outdir","$kaks_prog.results.txt"),
                    );

                $kaks_factory->alignment($dna_aln);
                my ($rc,$parser) = $kaks_factory->run();
                if ( $rc <= 0 ) { 
                    warn($kaks_factory->error_string,"\n");
                    exit;
                }
                my $result = $parser->next_result;
                my $MLmatrix = $result->get_MLmatrix();

                my @otus = $result->get_seqs();

                my @pos = map { 
                    my $c= 1;
                    foreach my $s ( @each) {
                        last if( $s->display_id eq $_->display_id );
                        $c++;
                    }
                    $c; 
                } @otus; 

                $out->_print
                    (
                     join("\t", qw(seq1 seq2 N S dN dS omega kappa lnL t cdna_pcntid)), "\n"
                    );
                for ( my $i = 0; $i < (scalar @otus -1) ; $i++) {
                    for ( my $j = $i+1; $j < (scalar @otus); $j++ ) {
                        my $sub_dna_aln = $dna_aln->select_noncont($pos[$i],$pos[$j]);
                        my $first_seq  = $otus[$i]->display_id;
                        my $second_seq = $otus[$j]->display_id;
                        my $N          = $MLmatrix->[$i]->[$j]->{'N'};
                        my $S          = $MLmatrix->[$i]->[$j]->{'S'};
                        my $dN         = $MLmatrix->[$i]->[$j]->{'dN'};
                        my $dS         = $MLmatrix->[$i]->[$j]->{'dS'};
                        my $omega      = $MLmatrix->[$i]->[$j]->{'omega'};
                        my $kappa      = $MLmatrix->[$i]->[$j]->{'kappa'};
                        my $lnL        = $MLmatrix->[$i]->[$j]->{'lnL'};
                        my $t          = $MLmatrix->[$i]->[$j]->{'t'};
                        my $percentid  = sprintf("%.2f",$sub_dna_aln->percentage_identity);
                        if (($dN == 0) || ($dS == 0)) {
                            $omega = 'nan';
                        }
                        if (($kappa == 999)) {
                            $kappa = 'nan';
                        }
                        $out->_print
                            (
                             join("\t",
                                  $first_seq ,
                                  $second_seq,
                                  $N        ,
                                  $S        ,
                                  $dN        ,
                                  $dS        ,
                                  $omega     ,
                                  $kappa     ,
                                  $lnL       ,
                                  $t         ,
                                  $percentid ,
                                 ), "\n"
                            );
                    }
                }
                $out->close();
            }
        }
    }
}

$logfile->close();

close DIR;
close SUBDIR;
close ALIGNDIR;

1;


