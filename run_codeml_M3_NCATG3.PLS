#!/usr/local/bin/perl

# run_codeml.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

run_codeml.PLS - DESCRIPTION 

=head1 SYNOPSIS

perl run_codeml.PLS \
-aln seqs.aln \
-tree species.tree \

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
use Bio::Tools::Run::Phylo::PAML::Codeml;
use Bio::AlignIO;
use Getopt::Long;

my ($aln_file, $tree_file,) = undef;
# my $outdir = undef;

GetOptions(
	   'a|aln:s' => \$aln_file,
	   't|tree:s' => \$tree_file,
          );
# 	   'o|outdir:s' => \$outdir,


BEGIN {$ENV{PAMLDIR} = '/home/contest/avb/9_opl/paml/paml3.14/src'; }

my $guessed_format = new Bio::Tools::GuessSeqFormat
    (-file=>"$aln_file"
    )->guess;

my $alignio = new Bio::AlignIO
    (
     -format=>$guessed_format,
     -file   => "$aln_file",
    );

my $treeio = new Bio::TreeIO
    (
     -format => 'newick', 
     -file => "$tree_file",
    );

my $tree = $treeio->next_tree;
my $aln = $alignio->next_aln;

my $codeml = new Bio::Tools::Run::Phylo::PAML::Codeml();
$codeml->alignment($aln);
$codeml->tree($tree);
$codeml->set_parameter("model","0");
#$codeml->set_parameter("model","0");
$codeml->set_parameter("noisy","9");
$codeml->set_parameter("verbose","0");
$codeml->set_parameter("runmode","0");
$codeml->set_parameter("seqtype","1");
$codeml->set_parameter("CodonFreq","2");
$codeml->set_parameter("clock","0");
$codeml->set_parameter("aaDist","0");
$codeml->set_parameter("aaRatefile","jones.dat");
$codeml->set_parameter("NSsites","3");
$codeml->set_parameter("icode","0");
$codeml->set_parameter("Mgene","0");
$codeml->set_parameter("fix_kappa","0");
$codeml->set_parameter("kappa","2");
$codeml->set_parameter("fix_omega","0");
$codeml->set_parameter("fix_alpha","1");
$codeml->set_parameter("alpha","0.");
$codeml->set_parameter("Malpha","0");
$codeml->set_parameter("ncatG","3");
$codeml->set_parameter("getSE","0");
$codeml->set_parameter("RateAncestor","1");
$codeml->set_parameter("Small_Diff",".5e-6");
$codeml->set_parameter("cleandata","1");
$codeml->set_parameter("fix_blength","1");
$codeml->set_parameter("method","0");
#$codeml->save_tempfiles(1);

print "Running codeml...\n";
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $start_time = sprintf ("%04d%02d%02d_%02d%02d%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec);
print " Started at $start_time\n";
my ($rc,$parser) = $codeml->run();
($sec,$min,$hour,$mday,$mon+1,$year,$wday,$yday,$isdst) = localtime(time);
my $end_time = sprintf ("%04d%02d%02d_%02d%02d%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec);
print "Finished at $end_time\n";
print "Parsing results...\n";
my $result = $parser->next_result;

my @trees = $result->get_trees;
#@trees = $result->get_rst_trees;

# Basic t w dN dS table for each branch
printf "anc_id;id;blen;branch;t;w;dN;dS;\n";
for my $t ( @trees ) {
    for my $node ( $t->get_nodes ) {	
        next unless $node->ancestor; # skip root node
        printf "%s;%s;%s;%s;%s;%s;%s;%s;\n",
            $node->ancestor->id,
            $node->id,
            $node->branch_length,
            $node->get_tag_values('branch'),
            $node->get_tag_values('t'),
            $node->get_tag_values('dN/dS'),
            $node->get_tag_values('dN'),
            $node->get_tag_values('dS');
    }
}

1;

################################################################################

# $codeml->set_parameter("","");
# 'noisy'   => [ 0..3,9],
# 'runmode' => [ -2, 0..5], 
# 'seqtype' => [ 1..3], # 1:codons, 2:AAs, 3:codons->AAs
# 'CodonFreq' => [ 2, 0,1,3], # 0:1/61 each, 1:F1X4, 
# 'aaDist'  => [ 0,'+','-', 1..6], 
# 'aaRatefile' => 'wag.dat', 
# 'model'    => [0..2,7], 
# 'NSsites'  => [0..13],
# 'icode'    => [ 0..10], 
# 'Mgene'    => [0,1], # 0:rates, 1:separate
# 'fix_kappa'=> [0,1], # 0:estimate kappa, 1:fix kappa
# 'kappa'    => '2',   # initial or fixed kappa
# 'fix_omega'=> [0,1], # 0: estimate omega, 1: fix omega
# 'omega'    => '0.4', # initial or fixed omega for 
# 'fix_alpha'=> [1,0], # 0: estimate gamma shape param
# 'alpha'    => '0', # initial of fixed alpha
# 'Malpha'   => [0,1], # different alphas for genes
# 'ncatG'    => [1..10], # number of categories in 
# 'clock'    => [0..3],
# 'getSE'    => [0,1], 
# 'RateAncestor' => [1,0,2], 
# 'Small_Diff'    => '.5e-6',
# 'cleandata'     => [0,1], 
# 'ndata'         => 1,
# 'method'        => [0,1], 
# 'fix_blength'   => [0,-1,2],

1;

