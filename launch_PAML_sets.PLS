#!/usr/local/bin/perl
# launch_PAML_sets.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

launch_PAML_sets.PLS - DESCRIPTION 

=head1 SYNOPSIS

Give standard usage here

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
use Bio::Tools::Phylo::PAML; # PAML parser
use Bio::Tools::Run::Phylo::PAML::Codeml;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Tools::GuessSeqFormat;
use Bio::Tools::Run::AnalysisFactory::Pise;
use Bio::Factory::EMBOSS;
# use Bio::Tools::CodonTable;

BEGIN {$ENV{PAMLDIR} = '/home/avb/9_opl/paml/paml3.14/src'; }
BEGIN {$ENV{EMBOSSDIR} = '/home/avb/bioperl-runnables/EMBOSS-2.9.0'; }

my ($factory, $application, $result, $emboss_present,
    $nucl_file, $aa_aln, $tranalignoutfile, 
    $guessed_format, $aligndir) = "";
my ($dir, $outdir, $subdir, $file, $input, $seq, 
    @seq_array, $seq_array_ref, $aln, $aa_seq,
    $codontable, $tableid, @params, $factory, $out) = "";

GetOptions(
 	   'table|tableid:s'=> \$tableid,
	   'dir:s' => \$dir,
	   'a|aligndir:s' => \$aligndir,
	   'o|outdir:s' => \$outdir,
          );

# get an EMBOSS application object from the EMBOSS factory
$factory = new Bio::Factory::EMBOSS;
$application = $factory->program('tranalign');

$emboss_present = $factory->executable();

unless ($emboss_present) {
    warn "emboss not found.\n";
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

            # running the application tranalign
            $tranalignoutfile = "$dir/$subdir/$outdir/$file.nucl";
            $tranalign->run
                (
                 {
                  '-asequence' => $nucl_file,
                  '-seqall'    => $aa_aln,
                  '-outseq'    => $tranalignoutfile,
                 }
                );

            # get the nucl alignment and run paml on it
            #  $alnin = new Bio::AlignIO
            #      (
            #       -format => $guessed_format,
            #       -file   => $tranalignoutfile
            #      );
            #  while ( my $aln = $alnin->next_aln ) {
            #      # process the alignment -- these will be Bio::SimpleAlign objects
            #  }

            $guessed_format = new Bio::Tools::GuessSeqFormat
                (
                 -file=>Bio::Root::IO->catfile("$dir","$subdir","$file")
                )->guess;
            $input = Bio::SeqIO->new
                (
                 -file=>Bio::Root::IO->catfile("$dir","$subdir","$file"),
                 -format=>$guessed_format,
                );
            while($seq = $input->next_seq()){
                $aa_seq = $seq->translate(undef, undef, undef, undef, $tableid);
                push (@seq_array, $aa_seq);
            }
            $seq_array_ref = \@seq_array;
            $aln = $factory->align($seq_array_ref);
            mkdir "$dir/$subdir/$outdir" or die "couldnt create directory: $!";
            $out = Bio::AlignIO->new(-file=>Bio::Root::IO->catfile(">","$dir","$subdir","$outdir","$file.phylip"),
                                     -format => 'phylip');
            $out->write_aln($aln);
            @seq_array = ();
        }
    }
}



my $codeml = new Bio::Tools::Run::Phylo::PAML::Codeml(-verbose => $verbose);
unless ($codeml->executable) {
  warn("PAML not is installed. skipping tests $Test::ntest to $NUMTESTS\n");
  exit(0) ;
}

my $in = new Bio::AlignIO(-format => 'phylip',
			  -file   => Bio::Root::IO->catfile(qw(t data 
							       gf-s85.phylip)));

my $aln = $in->next_aln;
$codeml->alignment($aln);
my ($rc,$results) = $codeml->run();

ok($rc,1);
if( ! defined $results ) { 
    exit(0);
}

my $result = $results->next_result;
my $MLmatrix = $result->get_MLmatrix;

my ($vnum) = ($result->version =~ /(\d+\.\d+)/);
# PAML 2.12 results
if( $vnum == 3.12 ) {
    ok($MLmatrix->[0]->[1]->{'dN'}, 0.0693);
    ok($MLmatrix->[0]->[1]->{'dS'},1.1459);
    ok($MLmatrix->[0]->[1]->{'omega'}, 0.0605);
    ok($MLmatrix->[0]->[1]->{'S'}, 273.5);
    ok($MLmatrix->[0]->[1]->{'N'}, 728.5);
    ok($MLmatrix->[0]->[1]->{'t'}, 1.0895);
} elsif( $vnum >= 3.13 ) {
# PAML 2.13 results
    ok($MLmatrix->[0]->[1]->{'dN'}, 0.0713);
    ok($MLmatrix->[0]->[1]->{'dS'},1.2462);
    ok(sprintf("%.4f",$MLmatrix->[0]->[1]->{'omega'}), 0.0572);
    ok($MLmatrix->[0]->[1]->{'S'}, 278.8);
    ok($MLmatrix->[0]->[1]->{'N'}, 723.2);
    ok(sprintf("%.4f",$MLmatrix->[0]->[1]->{'t'}), 1.1946);
} else { 
    for( 1..6) { 
	skip("Can't test the result output, don't know about PAML version ".$result->version,1);
    }
}

ok($codeml->error_string !~ /Error/); # we don't expect any errors;

my $yn00 = Bio::Tools::Run::Phylo::PAML::Yn00->new();
$yn00->alignment($aln);
($rc,$results) = $yn00->run();
ok($rc,1);
if( ! defined $results ) { 
    exit(0);
}
$result = $results->next_result;
$MLmatrix = $result->get_MLmatrix;

ok($MLmatrix->[0]->[1]->{'dN'}, 0.0846);
ok($MLmatrix->[0]->[1]->{'dS'}, 1.0926);
ok($MLmatrix->[0]->[1]->{'omega'}, 0.0774);
ok($MLmatrix->[0]->[1]->{'S'}, 278.4);
ok($MLmatrix->[0]->[1]->{'N'}, 723.6);
ok($MLmatrix->[0]->[1]->{'t'}, 1.0941);

ok($yn00->error_string !~ /Error/); # we don't expect any errors;


$codeml = new Bio::Tools::Run::Phylo::PAML::Codeml
    (-params => { 'alpha' => 1.53 },
     -verbose => $verbose);

ok($codeml);

1;


