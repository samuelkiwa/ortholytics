#!/usr/local/bin/perl
# add_GN_phylip_header.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

add_GN_phylip_header.PLS - DESCRIPTION 

=head1 SYNOPSIS

perl add_GN_phylip_header.PLS -i alignment_file.phy -h GN_header_file.txt

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

my ($inputfile,$headerfile,$yesnoopt);

GetOptions(
	   'i|input|inputfile:s' => \$inputfile,
           'h|header|headerfile:s' => \$headerfile,
          );

open HEADERFILE, "$headerfile" or die "$headerfile: $!";
open INPUTFILE, "$inputfile" or die "$inputfile: $!";
my $outfile = $inputfile;
$outfile =~ s|\.phy||g;
open OUTFILE, ">$inputfile.mgene.phy" or die "$inputfile.mgene.phy: $!";

my $lineheader;
while (<HEADERFILE>) {
    next unless ($_ =~ /G\d+/);
    $lineheader = $_;
}
close(HEADERFILE);


my $counter = 1;
while (<INPUTFILE>) {
    if (1 != $counter) {
        print OUTFILE $_;
    } else {
        my $g2line = $_;
        chomp $g2line;
        $g2line .= " G\n";
        $g2line .= "$lineheader";
        print OUTFILE $g2line;
    }
    $counter++;
}

close(INPUTFILE);
close(OUTFILE);


1;
