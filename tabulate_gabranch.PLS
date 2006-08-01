#!/usr/local/bin/perl
# tabulate_gabranch.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

tabulate_gabranch.PLS - DESCRIPTION 

=head1 SYNOPSIS

perl \
/home/avilella/ortholytics/tabulate_gabranch.PLS \
tabulate_gabranch.PLS \
-d \
/home/avilella/wallace/eukarya/drosophila/hc_ga_genes_wDgri_rlong

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

my ($dir,$yesnoopt);

GetOptions(
	   'd|dir:s' => \$dir,
           'yesnoopt' => \$yesnoopt,
          );

my $subdir;
my $file;
my $gabranch;
my $found_rates;
my $header = "seqfile,jobid,branchid,omega,ratevar,q025,q995,pbigw";
opendir DIR, $dir or die "couldnt find $dir:$!";
while (defined($subdir = readdir(DIR))) {
    next if ($subdir =~ /\./);
    next if ($subdir =~ /\.\./);
    opendir SUBDIR, "$dir/$subdir" or die "couldnt find subdir $subdir:$!";
    my $seqfile = `head -n 3 $dir/$subdir/$subdir.bf | tail -n 1`;
    $seqfile =~ s/.+\/(.+)\n/$1/;
    $found_rates = 0;
    while (defined($file = readdir(SUBDIR))) {
        #res.ga.s2.c2.r3.txt.samples.res_branch_equal.csv
        #res.ga.s2.c2.r3.txt.samples.res_branch_rates.csv
        if ($file =~ /(\S+)\.res_branch_equal.csv$/) {
#             open EQUAL,"$dir/$subdir/$file" or die "$!";
#             while (<EQUAL>) {
#                 1;
#                 $gabranch->{$subdir}{equal} = 1;
#             }
            1;
        } elsif ($file =~ /(\S+)\.res_branch_rates.csv$/) {
            open RATES,"$dir/$subdir/$file" or die "$!";
            while (<RATES>) {
                next if ($_ =~ /Branch\ Name/);
                #my ($branchid,$omega,$ratevar,$q025,$q995,$pbigw) = split(",", $_);
                $gabranch->{$subdir}{rates} .= "$seqfile,$subdir,". $_;
            }
            $found_rates = 1;
        }
    }
    unless ($found_rates) {
        print STDERR "Error reading rates file ", $subdir, "\n";
        1;
    }
}

print $header, "\n";
foreach my $subdir (sort keys %{$gabranch}) {
    print $gabranch->{$subdir}{rates};
}

1;
