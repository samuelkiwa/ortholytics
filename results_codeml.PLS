#!/usr/local/bin/perl

# results_codeml.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

results_codeml.PLS - DESCRIPTION 

=head1 SYNOPSIS

perl results_codeml.PLS
-mlcfile /my/paml/results/dir/mlc.txt
-dir /my/paml/results/dir

=head1 DESCRIPTION

results_codeml.PLS will load the results of a codeml mlc file and
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
use Bio::Tools::Phylo::PAML::ModelResult;
use Bio::SeqIO;

my ($mlcfile, $dir, $tag, $onlycodonfreqs,$fix);
$fix = 0;
$tag = 'notag';

GetOptions(
 	   'i|mlcfile|mlc|m:s'=> \$mlcfile,
	   'dir|d:s' => \$dir,
	   't|tag:s' => \$tag,
	   'fix' => \$fix,
	   'only|codon|freq' => \$onlycodonfreqs,
          );
my $parser = new Bio::Tools::Phylo::PAML
  (
   -file => "$mlcfile",
   -dir => "$dir",
  );
print STDERR "Parsing result $mlcfile...\n";

my $result = $parser->next_result();
my $tree = $result->next_tree;
my $seqfile = $result->seqfile();
$seqfile =~ s/.+\/(\S+)/$1/;
$seqfile =~ s/\.fasta|\.seqt\.phy//;

if ($onlycodonfreqs) {
    my @codon_freqs = $result->get_CodonFreqs();

    foreach my $firstbase (@codon_freqs) {
        foreach my $element (@$firstbase) {
            print "  $element";
        }
        print "\n";
    }
    exit 0;
}

# print "$tag;ancestorid..id;branch;t;S;N;dN/dS;dN;dS;S*dS;N*dN;\n";

# This will likely work for M0 vs FR runs
my @entries;
print STDERR "Entries will look like:\n";
print STDERR "tag,seqfile,branchid,lnL,bl,branchid2,t,S,N,omega,dN,dS,SdS,NdN\n";
push @entries,"tag,seqfile,branchid,lnL,bl,branchid2,t,S,N,omega,dN,dS,SdS,NdN\n";
eval {
    for my $node ( $tree->get_nodes ) {
        next unless $node->ancestor; # skip root node
        my $lnL = $tree->score;
        my $a = $node->ancestor->id;
        my $id = $node->id;
        my $namednode;
        if (length($id) < 2) {
            my $oldid = $id;
            $id = "Node";
            $id .= sprintf("%02d",$oldid);
        } elsif (length($a) < 2) {
            my $oldid = $a;
            $a = "Node";
            $a .= sprintf("%02d",$oldid);
            $namednode = $a;
        }
        my $bl = $node->branch_length;
        my $d = $node->get_tag_values('branch');
        my $e = $node->get_tag_values('t');
        my $f = $node->get_tag_values('S');
        my $g = $node->get_tag_values('N');
        my $h = $node->get_tag_values('dN/dS');
        $h = 5 if (($h > 5) && $fix);
        my $i = $node->get_tag_values('dN');
        my $j = $node->get_tag_values('dS');
        my $k = $node->get_tag_values('S*dS');
        my $l = $node->get_tag_values('N*dN');
        push @entries,"$tag,$seqfile,$id..$a,$lnL,$bl,$d,$e,$f,$g,$h,$i,$j,$k,$l\n",
    }
    foreach my $entry (@entries) {
        print "$entry";
    }
};

# This will likely work for M7 vs M8 runs
my @ns = $result->get_NSSite_results;
foreach my $ns (sort @ns) {
    my ($header,$values);
    my $time_used =         $ns->time_used;
    my $model_num =         $ns->model_num;
    my $model_description = $ns->model_description;
    my $kappa =             $ns->kappa;
    my $likelihood =        $ns->likelihood;
    my $tree =              $ns->next_tree;
#     print "time_used;        $time_used\n";
#     print "model_desc;       $model_description\n";
    $header .= "tag,";
    $values .= "$tag,";
    $header .= "seqfile,";
    $values .= "$seqfile,";
    $header .= "model_n,";
    $values .= "$model_num,";
    $header .= "kappa,";
    $values .= "$kappa,";
    $header .= "lnL,";
    $values .= "$likelihood,";
    my $num_site_classes;
    if ($num_site_classes = $ns->num_site_classes) {
        $header .= "nsite_classes,";
        $values .= "$num_site_classes,";
    }
    my %dnds_site_classes;
    if (%dnds_site_classes = %{$ns->dnds_site_classes}) {
        foreach my $key (sort keys %dnds_site_classes) {
            my $count = sprintf("%02d",1);
            foreach my $value (@{$dnds_site_classes{$key}}) {
                $header .= "$key"."i$count".",";
                $values .= "$value,";
                $count++; $count = sprintf("%02d",$count);
            }
            if (11 == $count) {
                $header .= "$key"."i$count".",";
                $values .= "0.00000,";
            }
        }
    }
    my $globalomega;
    for my $node ( $tree->get_nodes ) {
        $globalomega = $node->get_tag_values('dN/dS');
        next unless $node->ancestor; # skip root node
#         printf "%s;%s;%s;%s;%s;%s;%s;%s;\n",
#             $node->ancestor->id,
#             $node->id,
#             $node->branch_length,
#             $node->get_tag_values('branch'),
#             $node->get_tag_values('t'),
#             $node->get_tag_values('dN/dS'),
#             $node->get_tag_values('dN'),
#             $node->get_tag_values('dS');
    }
    $header .= "omega,";
    $values .= "$globalomega,";
    my %shape = %{$ns->shape_params};

    # Parameters in beta&w>1:
    #   p0=  1.00000  p=  0.07642 q=  0.85550
    #  (p1=  0.00000) w=  1.00000

    # Parameters in beta:
    #  p=  0.07642  q=  0.85549

    $header .= "betap,";
    $values .= "$shape{p},";
    $header .= "betaq,";
    $values .= "$shape{q},";

    $header .= "betap0,";
    if (defined $shape{p0}) {
        $values .= "$shape{p0},";
    } else {
        $values .= "1.00000,";
    }
    $header .= "betap1,";
    if (defined $shape{p1}) {
        $values .= "$shape{p1},";
    } else {
        $values .= "0.00000,";
    }
    $header .= "betaw";
    if (defined $shape{w}) {
        $values .= "$shape{w}";
    } else {
        $values .= "1.00000";
    }
    print "$header\n";
    print "$values\n";
}

1;
