#!/usr/local/bin/perl
# cross_go_ginput.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

cross_go_ginput.PLS - DESCRIPTION 

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
use Getopt::Long;

my ($listfile,$ginputfile,$yesnoopt);

GetOptions(
	   'l|listfile:s' => \$listfile,
	   'g|ginputfile:s' => \$ginputfile,
           'yesnoopt' => \$yesnoopt,
          );

open LISTFILE,"$listfile" or die "$!";
my %listentries;
my @genes;
while (<LISTFILE>) {
    if (/(\S+)\:(\S+)/) {
        @genes = split /\,/,$2;
        $listentries{$1} = @genes;
    } else {
        print "parsing error\n";
    }
}
open GINPUTFILE,"$ginputfile" or die "$!";
while (<GINPUTFILE>) {
    if (//) {

    } elsif (/branchid/) {
        # tag,seqfile,branchid,lnL,bl,branchid2,t,S,N,omega,dN,dS,SdS,NdN
    } else {
        print "parsing error\n";
    }
}

close LISTFILE;
close GINPUTFILE;

1;
