#!/usr/local/bin/perl

# codeml_results.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

codeml_results.PLS - DESCRIPTION 

=head1 SYNOPSIS

perl codeml_results.PLS
-mlcfile /my/paml/results/dir/mlc.txt
-dir /my/paml/results/dir

=head1 DESCRIPTION

codeml_results.PLS will load the results of a PAML mlc file and
directory (where rst and other companion files are found), and print
out the results in a spreadsheet-like format.

=head1 AUTHOR - Albert Vilella

Email 

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut


# Let the code begin...

use strict;
use Getopt::Long;
use Bio::Tools::Phylo::PAML;
use Bio::SeqIO;

my ($mlcfile, $dir);

GetOptions(
 	   'mlcfile|mlc|m:s'=> \$mlcfile,
	   'dir|d:s' => \$dir,
          );
my $parser = new Bio::Tools::Phylo::PAML
  (
   -file => "$mlcfile",
   -dir => "$dir",
  );
print "Parsing result...\n";

my $result = $parser->next_result();
my $tree = $result->next_tree;

# print "ancestorid..id;branch;t;S;N;dN/dS;dN;dS;S*dS;N*dN;\n";
for my $node ( $tree->get_nodes ) {
    next unless $node->ancestor; # skip root node
    my $a = $node->ancestor->id;
    my $b = $node->id;
    my $c = $node->branch_length;
    my $d = $node->get_tag_values('branch');
    my $e = $node->get_tag_values('t');
    my $f = $node->get_tag_values('S');
    my $g = $node->get_tag_values('N');
    my $h = $node->get_tag_values('dN/dS');
    my $i = $node->get_tag_values('dN');
    my $j = $node->get_tag_values('dS');
    my $k = $node->get_tag_values('S*dS');
    my $l = $node->get_tag_values('N*dN');
    print "$a..$b;$c;$d;$e;$f;$g;$h;$i;$j;$k;$l;\n";
}

my @ns = $result->get_NSSite_results;
foreach my $ns (sort @ns) {
    my $time_used =         $ns->time_used;
    my $model_num =         $ns->model_num;
    my $model_description = $ns->model_description;
    my $kappa =             $ns->kappa;
    my $likelihood =        $ns->likelihood;
    my $tree =              $ns->next_tree;
    print "time_used;        $time_used\n";
    print "model_n;          $model_num\n";
    print "model_desc;       $model_description\n";
    print "kappa;            $kappa\n";
    my $num_site_classes;
    if ($num_site_classes = $ns->num_site_classes) {
        print "nsite_classes; $num_site_classes\n";
    }
    my %dnds_site_classes;
    if (%dnds_site_classes = %{$ns->dnds_site_classes}) {
        foreach my $key (keys %dnds_site_classes) {
            print "$key;";
            foreach my $value (@{$dnds_site_classes{$key}}) {
                print "$value;";
            }
            print "\n";
        }
    }
    print "likelihood;$likelihood\n";
    for my $node ( $tree->get_nodes ) {
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

# my @seqs         = $result->get_seqs;
# my %input_params = $result->get_input_parameters;
# my @basfreq      = $result->get_codon_pos_basefreq;
# my $MLmatrix     = $result->get_MLmatrix; # get MaxLikelihood Matrix
# my $NGmatrix     = $result->get_NGmatrix; # get Nei-Gojoburi Matrix

1;;
