#!/usr/local/bin/perl
use Getopt::Long;

my ($inputfile,$outfile,$yesnoopt,$replicates,$seqgendir);

$seqgendir = "/home/avb/9_opl/seqgen/Seq-Gen.v1.3.2/source";
$pamldir = "/home/avb/9_opl/paml/_paml3.14b/paml3.14/src";

GetOptions(
	   't|tree|treefile:s' => \$inputfile,
	   'o|out|outfile:s' => \$outfile,
	   'r|repl:s' => \$replicates,
           'yesnoopt' => \$yesnoopt,
          );

while ($replicates) {
    $replicates = sprintf("%03d",$replicates);
    my $gamma_shape = int(rand(9)+2);
    my $kappa = sprintf("%.3f", rand(10));
    my $call = "$seqgendir"."/"."seq-gen -mHKY -g$gamma_shape -t$kappa -l1000 < $inputfile > $outfile.$replicates.fasta";
    system($call);
    print STDERR "# $call\n";
    $replicates--;
}



1;

