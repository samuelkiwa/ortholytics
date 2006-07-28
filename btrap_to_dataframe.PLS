#!/usr/local/bin/perl
# btrap_to_dataframe.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

btrap_to_dataframe.PLS - DESCRIPTION 

=head1 SYNOPSIS

perl btrap_to_dataframe.PLS \
-i \
/my/boot/file/res.Vest.1000.pboot.csv

=head1 DESCRIPTION

This script will create a r-project friendly csv file (dataframe) from
a typical hyphy free-ratios run with bootstrap.

=head1 AUTHOR - Albert Vilella

Email 

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut


# Let the code begin...

use strict;
use Getopt::Long;

my ($inputfile,$tag,$yesnoopt);

$tag = "tag";
GetOptions(
	   'i|input|inputfile:s' => \$inputfile,
	   'tag:s' => \$tag,
           'yesnoopt' => \$yesnoopt,
          );
my @tags = split /\:/,$tag;

open INFILE, "$inputfile" or die "cannot open $inputfile: $!";
my @columns;
my %entries;
my @dataframe;
my $csv_string = "branchid,branchtype,iter,lnL,t,omega,dN,dS";
my $counter = 1;
foreach my $mytag (@tags) {
    $csv_string .= ",";
    $csv_string .= "tag$counter";
    $counter++;
}

print STDERR "Entries will look like:\n";
print STDERR "$csv_string\n";
push @dataframe, $csv_string;

# Example header
# Iteration,Ln-likelihood,AC,givenTree.Node10.nonSynRate,givenTree.Node10.synRate,givenTree.P_cc9902.synRate,givenTree.P_cc9902.nonSynRate,givenTree.Pm_m9312.synRate,givenTree.P_cc9605.synRate,givenTree.Pm_m9312.nonSynRate,givenTree.Node9.nonSynRate,givenTree.Node8.nonSynRate,givenTree.Node8.synRate,givenTree.Pm_m9313.synRate,givenTree.Node9.synRate,givenTree.Pm_m9313.nonSynRate,givenTree.P_cc9605.nonSynRate,givenTree.Pm_natl2.nonSynRate,givenTree.P_wh8102.nonSynRate,givenTree.P_wh8102.synRate,givenTree.Node4.synRate,givenTree.Node2.nonSynRate,givenTree.Pm_cmed4.nonSynRate,givenTree.Pm_cmed4.synRate,givenTree.Pm_csarg.nonSynRate,givenTree.Pm_csarg.synRate,givenTree.Node2.synRate,givenTree.Node4.nonSynRate,givenTree.Pm_natl2.synRate,AT,CG,CT,GT,TotalLength(Pm_cmed4),SynLength(Pm_cmed4),NonSynLength(Pm_cmed4),TotalLength(Pm_natl2),SynLength(Pm_natl2),NonSynLength(Pm_natl2),TotalLength(P_wh8102),SynLength(P_wh8102),NonSynLength(P_wh8102),TotalLength(Node4),SynLength(Node4),NonSynLength(Node4),TotalLength(Node2),SynLength(Node2),NonSynLength(Node2),TotalLength(Pm_csarg),SynLength(Pm_csarg),NonSynLength(Pm_csarg),TotalLength(P_cc9605),SynLength(P_cc9605),NonSynLength(P_cc9605),TotalLength(Pm_m9312),SynLength(Pm_m9312),NonSynLength(Pm_m9312),TotalLength(Node10),SynLength(Node10),NonSynLength(Node10),TotalLength(P_cc9902),SynLength(P_cc9902),NonSynLength(P_cc9902),TotalLength(Node9),SynLength(Node9),NonSynLength(Node9),TotalLength(Pm_m9313),SynLength(Pm_m9313),NonSynLength(Pm_m9313),TotalLength(Node8),SynLength(Node8),NonSynLength(Node8)

while (<INFILE>) {
    # Header -- first line
    if ($_ =~ /^Iteration/) {
        @columns = split /\,/,$_;
        my $counter = 0;
        foreach my $colname (@columns) {
            if ($colname =~ /Length/ ||
                $colname =~ /likelihood/ ||
                $colname =~ /Iteration/) {
                $colname =~ s/\(/\./;
                $colname =~ s/\)//;
                $colname =~ s/Ln-likelihood/lnL/;
                $entries{$colname} = $counter;
            }
            $counter++;
        }
    # Data line
    } else {
        @columns = split /\,/,$_;
        foreach my $entry (sort keys %entries) {
            if ($entry =~ /TotalLength/) {
                $entry =~ /Length\.(.+)/;
                my $branch_id = $1;
                my $csv_string = "$branch_id";
                my $branch_type;
                if ($branch_id =~ /Node/) {
                    $branch_type = "internal";
                } else {
                    $branch_type = "external";
                }
                $csv_string .= ",";
                my $dN = "$columns[$entries{\"NonSynLength.$branch_id\"}]";
                my $dS = "$columns[$entries{\"SynLength.$branch_id\"}]";
                $dN = sprintf("%.7f",$dN);
                $dS = sprintf("%.7f",$dS);
                my $omega;
                if ($dS != 0) {
                    $omega = sprintf("%.7f",$dN/$dS);
                } else {
                    $omega = sprintf("%.7f",3);
                }
                $csv_string .= "$branch_type";
                $csv_string .= ",";
                $csv_string .= "$columns[$entries{Iteration}]";
                $csv_string .= ",";
                $csv_string .= "$columns[$entries{\"lnL\"}]";
                $csv_string .= ",";
                $csv_string .= "$columns[$entries{\"TotalLength.$branch_id\"}]";
                $csv_string .= ",";
                $csv_string .= "$omega";
                $csv_string .= ",";
                $csv_string .= "$dN";
                $csv_string .= ",";
                $csv_string .= "$dS";
                foreach my $mytag (@tags) {
                    $csv_string .= ",";
                    $csv_string .= "$mytag";
                }
                push @dataframe, $csv_string;
            }
        }
    }
}

my $outputfile = "$inputfile.dataframe.csv";
open OUTFILE, ">$outputfile" or die "cannot open $outputfile: $!";
foreach my $entry (@dataframe) {
    print OUTFILE "$entry\n";
}
close(OUTFILE);
1;
