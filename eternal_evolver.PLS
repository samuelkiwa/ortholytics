# eternal_evolver.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

eternal_evolver.PLS - DESCRIPTION 

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
use Bio::Tools::Run::Phylo::PAML::Evolver;
use Bio::Tools::Run::Phylo::PAML::Codeml;
use Getopt::Long;
use File::Path;

my ($workdir, $exedir, $seqfile, $treefile, $branchModel, $nssites,
    $codonfreq, $class, $mgene, $sendmail, $nomail, $zip, $verbose,
   $rateancestor, $ncatG) =
    undef;

my ($runid,$program) = undef;
$branchModel = 0; $nssites = 0; $codonfreq = 5;
$runid = 1;
$program = "codeml";
$class = "research_long";
$mgene = 0;
## this two enable further evolver stuff
$verbose = 1;
$rateancestor = 1;
$ncatG = 3;

GetOptions(
	   'w|work|workdir:s'  => \$workdir,
	   'e|exe|exedir:s'    => \$exedir,
	   's|seq|seqfile:s'   => \$seqfile,
	   't|tree|treefile:s' => \$treefile,
           'b|branch|branch_model|branchModel|model:s' => \$branchModel,
	   'n|ns|nssites:s'    => \$nssites,
	   'c|cod|codonfreq:s' => \$codonfreq,
	   'class|queue:s'    => \$class,
	   'mgene:s'    => \$mgene,
	   'ncatg:s'    => \$ncatG,
	   'zip|zipped'    => \$zip,
           'verbose:s' => \$verbose,
           'a|rateancestor:s' => \$rateancestor,
	   'nomail'    => \$nomail,
	   'sendmail' => \$sendmail,
          );

BEGIN {$ENV{PAMLDIR} = '/home/avb/9_opl/paml/_paml3.14b/paml3.14/src'; };

my $analysis;
my ($dir, $batch, $loader);
$runid = sprintf("%04d", $runid);
$analysis = "$program"."$runid";
$dir = "$workdir/$analysis";
while (-d $dir) {
    $runid++;
    $runid = sprintf("%04d", $runid);
    $analysis = "$program"."$runid";
    $dir = "$workdir/$analysis";
}

unless (-d $dir) {
    File::Path::mkpath($dir);
#     mkdir $dir or die "couldnt create directory: $!";
}


my $alignio = new Bio::AlignIO
    (
     -format => "phylip",
     -file   => "$seqfile",
    );

my $treeio = new Bio::TreeIO
    (
     -format => 'newick', 
     -file => "$treefile",
    );

my $tree = $treeio->next_tree;
my $aln = $alignio->next_aln;

my $codeml = new Bio::Tools::Run::Phylo::PAML::Codeml();
# $codeml->executable("$exedir");
$codeml->alignment($aln);
$codeml->tree($tree);

$codeml->no_param_checks(1);
if ($branchModel =~ /m0/i) {$branchModel = "0";} 
elsif ($branchModel =~ /fr/i) {$branchModel = "1";}
$codeml->set_parameter("model","$branchModel");
$codeml->set_parameter("noisy","9");
$codeml->set_parameter("verbose","$verbose");
$codeml->set_parameter("runmode","0");
$codeml->set_parameter("seqtype","1");
$codeml->set_parameter("CodonFreq","$codonfreq");
$codeml->set_parameter("clock","0");
$codeml->set_parameter("aaDist","0");
$codeml->set_parameter("aaRatefile","jones.dat");
$codeml->set_parameter("NSsites","$nssites");
$codeml->set_parameter("icode","0");
$codeml->set_parameter("Mgene","$mgene");
$codeml->set_parameter("fix_kappa","0");
$codeml->set_parameter("kappa","2");
$codeml->set_parameter("fix_omega","0");
$codeml->set_parameter("fix_alpha","1");
$codeml->set_parameter("alpha","0.");
$codeml->set_parameter("Malpha","0");
$codeml->set_parameter("ncatG","$ncatG");
$codeml->set_parameter("getSE","0");
$codeml->set_parameter("RateAncestor","$rateancestor");
$codeml->set_parameter("Small_Diff",".5e-6");
$codeml->set_parameter("cleandata","1");
$codeml->set_parameter("fix_blength","1");
$codeml->set_parameter("method","0");
$codeml->save_tempfiles(1);

$codeml->tempdir($dir);
my ($rc,$results) = $codeml->run();

# Create an evolver preparation instance in the specified workdir
my $evolver = new Bio::Tools::Run::Phylo::PAML::Evolver();
unless (-d $workdir) {
    eval "require File::Path";
    if ($@) {
        print "File::Path not found. trying with mkdir\n";
        mkdir("$workdir");
    } else {
        File::Path::mkpath($workdir);
    }
};
$evolver->tempdir("$workdir");
$evolver->save_tempfiles(1);

my $inputfile = "$dir/mlc";

# With all the data coming from a previous PAML run
my $parser = new Bio::Tools::Phylo::PAML
  (
   -file => "$inputfile",
   -dir => "$dir",
  );
my $result = $parser->next_result();

# Codon frequencies as in the codeml run
my @codon_freqs = $result->get_CodonFreqs();
$evolver->set_CodonFreqs(\@codon_freqs);

# Check the length of the rst seqs. This length is most of the times
# different from the original data due to gaps in the MSA
my @rsts = $result->get_rst_seqs;
my $nuclsites = $rsts[0]->length;
$nuclsites = $nuclsites/3;

$evolver->set_parameter("nuclsites","$nuclsites");
$evolver->set_parameter("tree_length","1");

# Get the kappa parameter from an NSSite results This needs a double
# codeml run (e.g., nssites ="0 1") to be parsed
my @ns = $result->get_NSSite_results;
my $kappa;
eval {$kappa = $ns[0]->kappa;};
if ($@) {
    $evolver->set_parameter("kappa","2");
} else {
    $evolver->set_parameter("kappa","$kappa");
}

# Tree from NSSites
# eval {
#     $tree =  $ns[0]->next_tree;
# };
$evolver->tree($tree);

#FIXME: get a proper omega, for God sake
my $NGmatrix = $result->get_NGmatrix;
my $dummyomega = $NGmatrix->[0]->[1]->{omega};
$evolver->set_parameter("omega","$dummyomega");
#$evolver->set_parameter(omega,"3\n0.2\t0.3\t0.5\n0.5\t0.9\t3.2\n");

# 3            * number of site classes, followed by frequencies and omega's.
#   0.6    0.3   0.1 # Freqs
#   0.1    0.8   3.2 # Omegas

# All ready - prepare
# $evolver->executable("$exedir");
my $dummy = $evolver->prepare();
my $rc = $evolver->run();

1;;
