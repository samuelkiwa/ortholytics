# evolver_indel.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

evolver_indel.PLS - DESCRIPTION 

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

use Getopt::Long;

my ($treefile,$dir,$rescale,$fakerep,$debug);

GetOptions(
	   'i|tree|treefile:s' => \$treefile,
	   'dir:s' => \$dir,
           'debug' => \$debug,
           'rescale' => \$rescale,
           'fakereplicates:s' => \$fakerep,
          );

my $fakereplicates = $fakerep || 1;

BEGIN {$ENV{EVOLVER_INDELDIR} = '/home/avb/9_opl/evolver_indel/evolver_indel'; };
#BEGIN {$ENV{PAMLDIR} = '/home/avilella/9_opl/paml/_paml3.14b/paml3.14/src'; };

my $treeio = new Bio::TreeIO
    (
     -format => 'newick', 
     -file => "$treefile",
    );

my $tree = $treeio->next_tree;

my $total_branch_length = 0;
if ($rescale) {
    foreach my $node ($tree->get_nodes) {
        $total_branch_length += $node->branch_length;
    }
    foreach my $node ($tree->get_nodes) {
        my $branch_length = $node->branch_length;
        next unless (defined($branch_length));
        $node->branch_length(sprintf("%.10f",$branch_length/$total_branch_length));
    }
    my $new_branch_length;
    foreach my $node ($tree->get_nodes) {
        $new_branch_length += $node->branch_length;
    }
}
if ($debug) {
    my $treeout = Bio::TreeIO->new(-format => 'newick'); #print to STDOUT instead
    $treeout->write_tree($tree);
}

foreach my $replicate (1 .. $fakereplicates) {
    my $repldir = "$dir/$replicate";
    unless (-d $repldir) {
        eval "require File::Path";
        if ($@) {
            print "File::Path not found. trying with mkdir\n";
            mkdir("$repldir");
        } else {
            File::Path::mkpath($repldir);
        }
    }
    my $evolver = new Bio::Tools::Run::Phylo::PAML::Evolver();
    $evolver->tempdir("$dir");
    $evolver->save_tempfiles(1);
    my $dummynuclsites = (int(rand(10000)))+1;
    $evolver->set_parameter("nuclsites","$dummynuclsites");
    $evolver->set_parameter("tree_length","1");
    $evolver->set_parameter("kappa","2");
    #FIXME - tree to sum 1
    $evolver->tree($tree);
    my $treeoutfile = Bio::TreeIO->new(-format => 'newick', -file => ">$dir/evolver_indel.input_tree.nh"); #print to STDOUT instead
    $treeoutfile->write_tree($tree);
    my $dummyomega = (int(rand(10))/10);
    $evolver->set_parameter("omega","$dummyomega");
    $evolver->tempdir("$repldir");
    $evolver->save_tempfiles(1);
    my $dummy = $evolver->prepare();
    my $rc = $evolver->run();
    print $evolver->error_string, "\n";
}

1;
