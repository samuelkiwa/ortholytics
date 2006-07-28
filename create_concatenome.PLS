# create_concatenome.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

perl create_concatenome.PLS -dir /my/dir -aligndir probcons

=head1 SYNOPSIS

This script will read the nucleotide MSA (dna.fasta) files from each
subdirectory (subdir, a bunch of them, all in maxdepth 1) in aligndir,
and concatenate them all, to form a big MSA (dna.fasta). All the
directory names start in an absolute path defined in dir option. It
will print in the log file the paml header necessary for the "Option
G2", i.e., with the codon lengths of the genes in the concatenome.

It will generate an aln stats file with, among others:

average_percentage_identity: This method implemented by Kevin Howe
                             calculates a figure that is designed to
                             be similar to the average pairwise
                             identity of the alignment (identical in
                             the absence of gaps), without having to
                             explicitly calculate pairwise identities
                             proposed by Richard Durbin.  Validated by
                             Ewan Birney ad Alex Bateman.

overall_percentage_identity: The percentage identity of the conserved
                             columns


=head1 DESCRIPTION

perl create_concatenome.PLS -dir /my/dir -aligndir probcons 

[-aa]         aminoacidic concatenome (instead of nucleotide concat)

[-th 50]      only a 50th of the genes will be concatenated

[-tag foo]    add a foo tag in the resulting name file

[-highlimit 95.0] put a high threshold to the overall percentage
                  identity of each aln to be included in the
                  concatenome. Ex: include only alignments with
                  identity of 95.0 or less.

[-lowlimit 20.0]  put a low threshold to the overall percentage
                  identity of each aln to be included in the
                  concatenome. Ex: include only alignments with
                  identity of 20.0 or more.

[-aa]         give aminoacidic concatenome instead of nucleotide

[-namelist] comma-separated list of the names of the sequences (in the
              same order as they are alphabetically in the input
              alignments)

[-genelist] a list of genes


=head1 AUTHOR - Albert Vilella

Email 

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut


# Let the code begin...

use strict;

use Getopt::Long;
use Bio::AlignIO;
use Bio::Align::DNAStatistics;

my ($aln_file, $aln_input, $aln, %seqs, $newAlignment, 
    $length, $g2_string, $out, $out_aln);

my ($dir, $subdir, $aligndir, $outdir_opt, 
    $logfile, $gnum_header_file, $gl_file, $stats_file, 
    $email, $highlimit_opt, $lowlimit_opt, $aminoconcat, 
    $th, $tag, $timetag, $namelist, $genelist,
    $dostats, $dokimura) = "";
$dostats = 0;
$dokimura = 0;
$th = 1;
my @opt_string = @ARGV;
GetOptions(
	   'd|dir:s' => \$dir,
	   'a|aligndir:s' => \$aligndir,
	   'o|outdir:s' => \$outdir_opt,
           'l|low|lowlimit:s' => \$lowlimit_opt,
           'h|high|highlimit:s' => \$highlimit_opt,
           'aa|amino|aminoconcat' => \$aminoconcat,
           't|tag:s' => \$tag,
           'th:s' => \$th,
           'n|namelist|names:s' => \$namelist,
           'g|genelist:s' => \$genelist,
           'timetag' => \$timetag,
           'dostats' => \$dostats,
           'dokimura' => \$dokimura,
	   'e|email:s' => \$email,
          );

my $highlimit = $highlimit_opt || 100;
my $lowlimit = $lowlimit_opt || 0;
my @names = split /\:/,$namelist if ($namelist);
my $outdir = $outdir_opt || $dir;
unless (-d $outdir) {
    eval "require File::Path";
    if ($@) {
        print "File::Path not found. trying with mkdir\n";
        mkdir("$outdir");
    } else {
        File::Path::mkpath($outdir);
    }
};

print "Concatenate alignments in a big concatenome...\n";
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $logfile_time = 
    sprintf ("%04d%02d%02d_%02d%02d%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec);

my $outfile_name = "concatenome";
if ($genelist) {
    my $genelist_id = $genelist;
    $genelist_id =~ /(\w+\d+)$/;
    $genelist_id = $1;
    $tag = $genelist_id unless($tag);
}
if ($aminoconcat)     {    $outfile_name .= ".aa"; }
if ($lowlimit_opt)    {    $outfile_name .= ".idl" . "$lowlimit"; }
if ($highlimit_opt)   {    $outfile_name .= ".idh" . "$highlimit"; }
if ($tag)             {    $outfile_name .= ".$tag"; }
if ($timetag)         {    $outfile_name .= ".$logfile_time"; }

$logfile = Bio::Root::IO->new
    (-file=>Bio::Root::IO->catfile
     (">","$outdir","logfile.$outfile_name.log"
     ),);

$gnum_header_file = Bio::Root::IO->new
    (-file=>Bio::Root::IO->catfile
     (">","$outdir","gnum_header.$outfile_name.log"
     ),);

$gl_file = Bio::Root::IO->new
    (-file=>Bio::Root::IO->catfile
     (">","$outdir","gene_lengths.$outfile_name.log"
     ),);

$stats_file = Bio::Root::IO->new
    (-file=>Bio::Root::IO->catfile
     (">","$outdir","aln_stats.$outfile_name.log"
     ),);

my $firstseq = 0; my $num_processed = 0;
$out = new Bio::SimpleAlign();
$out_aln = new Bio::AlignIO
    (
     -file=>Bio::Root::IO->catfile(">", "$outdir","$outfile_name.fasta"),
     -format=> 'fasta',
    );

printf ("Writing log to %s\n", $logfile->file);

$logfile->_print
    ("Started $logfile_time\n");

$logfile->_print
    ("With options: @opt_string\n");

my @genelist;
if ($genelist) {
    open(GENELIST,"$genelist") or die "couldnt open file $genelist: $!\n";
    while (<GENELIST>) {
        chomp;
        next if ($_ =~ /^#/);
        my @fields = split /\,/,$_;
        push @genelist,$fields[0];
    }
}

my $thconcatenome = 0;

if ($dokimura) {
    $stats_file->_print("num_processed;aln_file;no_sequences;aln_length;transitions;transversions;avg_k2p;percentage_identity;cnsv_cols_percentage_identity;\n");
} else {
    $stats_file->_print("numprocessed;alnfile;nosequences;alnlength;transitions;transversions;percentageidentity;cnsvcolspercentageidentity;\n");
}

opendir DIR, $dir or die "couldnt find $dir:$!";
while (defined($subdir = readdir(DIR))) {
    next if ($subdir =~ /\./);
    next if ($subdir =~ /\.\./);
    opendir ALIGNDIR, "$dir/$subdir/$aligndir" or die "couldnt find $dir/$subdir/$aligndir:$!";
    while (( defined($aln_file = readdir(ALIGNDIR)) )) {
        #Look for the alignments in the dir
        if ($aln_file =~ /(\w+)(\-\w+)?\.fasta\.afa(\.dna.fasta)?$/
            ||
           $aln_file =~ /mc\.paml/
            ||
           $aln_file =~ /ancestral\.seq/) {
            my $not_in_list = 0;
            if ($genelist) {
                $not_in_list = 1 unless ( (grep { /^$1$/i } @genelist) ) ;
            }
            next if ($not_in_list);
            my $guessed_format = new Bio::Tools::GuessSeqFormat
                (-file=>Bio::Root::IO->catfile("$dir","$subdir","$aligndir","$aln_file")
                )->guess;
            $aln_input = Bio::AlignIO->new
                (
                 -file=>Bio::Root::IO->catfile("$dir","$subdir","$aligndir","$aln_file"),
                 -format=>$guessed_format,
                );

            $aln = $aln_input->next_aln();

            # Run aln statistics
            my $stats = new Bio::Align::DNAStatistics;
            my $cnsv_cols_percentage_identity = 0;
            $cnsv_cols_percentage_identity = $aln->overall_percentage_identity if ($dostats);

            # first alignment is the base to which others are added
            my ($toaddSeq, $concatSeq, $id);

            if ($cnsv_cols_percentage_identity >= $lowlimit && $cnsv_cols_percentage_identity <= $highlimit) {
                $thconcatenome++;
                next unless (!($thconcatenome%$th));
                $num_processed++;
                $logfile->_print("no_$num_processed;$aln_file;$dir/$subdir/$aligndir/$aln_file;\n");
                $gl_file->_print("$num_processed;$dir/$subdir/$aligndir/$aln_file/;$subdir;");

                my $transitions   = 0;
                my $transversions = 0;
                unless ($aminoconcat) {
                     $transitions   = $stats->transitions($aln) if ($dostats);
                     $transversions = $stats->transversions($aln) if ($dostats);
                }
                my $percentage_identity = 0;
                $percentage_identity = $aln->percentage_identity if ($aminoconcat);
                my $average_percentage_identity = 0;
                $average_percentage_identity = $aln->average_percentage_identity if ($aminoconcat);
                my $no_sequences                = $aln->no_sequences;
                my $aln_length                  = $aln->length;

                if ($dostats && $dokimura) {
                    my $k2pmatrix = $stats->distance(-align => $aln, 
                                                     -method => 'Kimura');
                    my @distances = split /\s+/, $k2pmatrix->print_matrix;
                    my $avg_k2p = 0;
                    my $numc = 0;
                    foreach my $entry (@distances) {
                        if ($entry =~ /[0-9]*\.[0-9]*/) {
                            unless (0.0 >= $entry) {
                                $avg_k2p += $entry;
                                $numc++;
                            }
                        }
                    }
                    $avg_k2p = $avg_k2p/($numc);
                    my ($Nei_Gojobori,$D_s,$D_n,$D_s_var,$D_n_var) = undef;
                    eval {
                        $Nei_Gojobori = $stats->calc_average_KaKs($aln, 10);
                        $D_s = $Nei_Gojobori->{D_s};
                     $D_n = $Nei_Gojobori->{D_n};
                     $D_s_var = $Nei_Gojobori->{D_s_var};
                     $D_n_var = $Nei_Gojobori->{D_n_var};
                         };
                        $stats_file->_print("$num_processed;$aln_file;$no_sequences;$aln_length;$transitions;$transversions;$avg_k2p;$percentage_identity;$cnsv_cols_percentage_identity;\n");
                } elsif ($dostats && !$dokimura) {
                    $stats_file->_print("$num_processed;$aln_file;$no_sequences;$aln_length;$transitions;$transversions;$percentage_identity;$cnsv_cols_percentage_identity;\n");
                }
                if ($firstseq == 0) {
                    my $count=1;
                    $aln->sort_alphabetically;
                    my $numseq = $aln->no_sequences;
                    for ($count = 1; $count <= $numseq; $count++) {
                        $toaddSeq = $aln->get_seq_by_pos($count);
                        $seqs{$count} .= $toaddSeq->seq;
                        $length = $toaddSeq->length;
                        $logfile->_print("length $length;");
                        if (($count == 1) && !($aminoconcat)) {
                            my $codons_length = $length/3;
                            $g2_string .= "GNUM\t";
                            $g2_string .= "$codons_length ";
                            $gl_file->_print("$length;");
                        }
                    }
                    $firstseq = 1;

                } else {
                    my $count=1;
                    $aln->sort_alphabetically;
                    my $numseq = $aln->no_sequences;
                    for ($count = 1; $count <= $numseq; $count++) {
                        $toaddSeq = $aln->get_seq_by_pos($count);
                        $seqs{$count} .= $toaddSeq->seq;
                        $length = $toaddSeq->length;
                        $logfile->_print("length $length;");
                        my $accum_length = length $seqs{$count};
                        $logfile->_print("accum_length $accum_length;");
                        if (($count == 1)  && !($aminoconcat)) {
                            my $codons_length = $length/3;
                            $g2_string .= "$codons_length ";
                            $gl_file->_print("$length;");
                        }
                    }
                }
                $gl_file->_print("\n");
            }
        }
    }
}

$logfile->_print("Final write_aln");
foreach my $seq_key (keys %seqs) {
    my $id = sprintf("%03d", $seq_key);
    my $length = length $seqs{$seq_key};
    $logfile->_print(" final length $length\n");
    my $long_id = $names[($seq_key-1)] || "sp_$id";
    my $newSeq = new Bio::LocatableSeq(-seq => $seqs{$seq_key},
                                       -id  => "$long_id",
#                                        -start => 1,
#                                        -end   => $length
                                      );
    $out->add_seq($newSeq);
    $seqs{$seq_key} = "";
}

$g2_string =~ s/GNUM/G$num_processed/g;
unless ($aminoconcat) {
    $gnum_header_file->_print
        (
         "\n$g2_string\n"
        );
}

$out->set_displayname_flat(1);
$out_aln->write_aln($out);
$out_aln->close();
$logfile->close();
$gnum_header_file->close();
$stats_file->close();
$gl_file->close();

close DIR;
close SUBDIR;
close ALIGNDIR;

################################################################################
# Email results via Mail::Sendmail

my $start_time = $logfile_time;
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $end_time = sprintf ("%04d%02d%02d_%02d%02d%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec);

my $message = "Program finished.";
$message .= " \n\nStarted at $start_time\nFinished at $end_time\n\n";

if ($email) {
    eval "require Mail::Sendmail";
    if ($@) {
        print "Mail::Sendmail not found. Email wont be sent\n";
    } else {
        my %mail = ( To      => "$email",
                     From    => 'create_concatenome@localhost',
                     Subject => "Sendmail batch: job finished.",
                     Message => "$message"
                   );

        sendmail(%mail) or die $Mail::Sendmail::error;
        print "Mail sent. Log says:\n", $Mail::Sendmail::log;
        print "\n";
    }
}

1;
