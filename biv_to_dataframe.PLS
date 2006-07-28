#!/usr/local/bin/perl
# biv_to_dataframe.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

biv_to_dataframe.PLS - DESCRIPTION 

=head1 SYNOPSIS

perl biv_to_dataframe.PLS \
-i \
/my/biv/file/res.bivariate.6.summary.txt

=head1 DESCRIPTION

This script will create a r-project friendly csv file (dataframe) from
a typical hyphy bivariate run.

=head1 AUTHOR - Albert Vilella

Email 

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut


# Let the code begin...

use strict;
use Getopt::Long;

my ($inputfile,$tag,$quiet);

$tag = "notag";
GetOptions(
	   'i|input|inputfile:s' => \$inputfile,
	   'tag:s' => \$tag,
           'quiet' => \$quiet,
          );
my @tags = split /\:/,$tag;
my $counter = 1;
my $csv_string = "";
foreach my $mytag (@tags) {
    $csv_string .= "tag$counter";
    $csv_string .= ",";
    $counter++;
}
$tag =~ s/\:/\,/g;

open INFILE, "$inputfile" or die "cannot open $inputfile: $!";
my %entries;
my @dataframe;

$csv_string .= "class,lnL,params,aic,caic,dS,dN,omega,prob\n";
print STDERR "Entries will look like:\n" unless ($quiet);
print STDERR "$csv_string\n" unless ($quiet);
push @dataframe, $csv_string;

# Model fit summary

# Log likelihood:    -1184.61321
# Parameters    :             45
# AIC           :     2459.22643
# c-AIC         :     2540.40290
# Class 1
#         dS    =      1.179
#         dN    =      0.262
#         dN/dS =      0.222
#         Prob  =      0.271
# Class 2
#         dS    =      0.934
#         dN    =      0.014
#         dN/dS =      0.015
#         Prob  =      0.727
# Class 3
#         dS    =      0.552
#         dN    =      0.552
#         dN/dS =      1.000
#         Prob  =      0.001
# Class 4
#         dS    =      0.549
#         dN    =      0.549
#         dN/dS =      1.000
#         Prob  =      0.000
# Class 5
#         dS    =      0.881
#         dN    =      1.763
#         dN/dS =      2.000
#         Prob  =      0.000
# Class 6
#         dS    =      0.881
#         dN    =      2.644
#         dN/dS =      3.000
#         Prob  =      0.000

my $class;
my ($lnL,$params,$aic,$caic);
while (<INFILE>) {
    next if $_ =~ /^$/;
    next if $_ =~ /^Model fit summary/;
    if ($_ =~ /likelihood\:\s*(\S+)/) {
        $lnL = $1;
    } elsif ($_ =~ /Parameters\s*\:\s*(\S+)/) {
        $params = $1;
    } elsif ($_ =~ /^AIC\s*\:\s*(\S+)/) {
        $aic = $1;
    } elsif ($_ =~ /^c-AIC\s*\:\s*(\S+)/) {
        $caic = $1;
    } else {
        if ($_ =~ /Class\s*(\S+)/) {
            my $num = sprintf("%02d", $1);
            $class = "c" . "$num";
        } elsif ($_ =~ /\s*(\S+)\s*\=\s*(\S+)/) {
            $entries{$class}{$1} = $2;
        }
    }
}
foreach my $class (sort keys %entries) {
    $csv_string = "$tag,$class,$lnL,$params,$aic,$caic,$entries{$class}{'dS'},$entries{$class}{'dN'},$entries{$class}{'dN/dS'},$entries{$class}{'Prob'}\n";
    push @dataframe, $csv_string;
}

my $outputfile = "$inputfile.dataframe.csv";
open OUTFILE, ">$outputfile" or die "cannot open $outputfile: $!";
foreach my $entry (@dataframe) {
    print OUTFILE "$entry";
}
close(OUTFILE);
1;
