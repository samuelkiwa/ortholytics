# evolver_rsts.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

evolver_rsts.PLS - DESCRIPTION 

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
use Getopt::Long;

my ($inputfile,$dir,$outdir,$yesnoopt);

GetOptions(
	   'm|i|input|inputfile:s' => \$inputfile,
	   'd|dir:s' => \$dir,
	   'o|outdir:s' => \$outdir,
           'yesnoopt' => \$yesnoopt,
          );

# Create an evolver preparation instance in the specified outdir
my $evolver = new Bio::Tools::Run::Phylo::PAML::Evolver();
unless (-d $outdir) {
    eval "require File::Path";
    if ($@) {
        print "File::Path not found. trying with mkdir\n";
        mkdir("$outdir");
    } else {
        File::Path::mkpath($outdir);
    }
};
$evolver->tempdir("$outdir");
$evolver->save_tempfiles(1);

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
my $tree =  $ns[0]->next_tree;
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
my $dummy = $evolver->prepare();

# A simple prompt so that everything can be executed now
print STDERR "# Creating... ",$evolver->tempdir,"\n";
print STDERR 
    " cd ",$evolver->tempdir,
    " && ",
    "\$PAMLDIR/evolver 6 MCcodon.dat",
    " && ",
    "perl -i.old -ne 'print unless 1 .. 6' ancestral.seq",
    " && ",
    "cd -\n";

print STDERR 
    "# cd ",$evolver->tempdir,
    " && ",
    "\$PAMLDIR/evolverNSsites 6 MCcodon.dat",
    " && ",
    "perl -i.old -ne 'print unless 1 .. 6' ancestral.seq",
    " && ",
    "cd -\n";

1;;
