#!/usr/bin/env perl
#
# Copyright (c) 2004 Alexander Anderson. All rights reserved.
# This program is free software; you can redistribute it and/or
# modify it under the same terms as Perl itself.
#
# $Id: xls2csv.pl,v 1.7 2005/01/07 15:46:49 anderal Exp anderal $

use Getopt::Std;
use Spreadsheet::ParseExcel::Simple;
use Fatal qw(open);

use strict;

our($opt_s, $opt_v);
getopts('vs:');

my %indexes = map { $_ => 1 } split(/,/, $opt_s);

# Doing it one Excel file at a time because Spreadsheet::ParseExcel seems to
# have a memory leak when dealing with lots (>12) large (>1Mb) Excel files.
# (or maybe memory leak used to be in the script?)

my $filename = shift or die "Usage: xls2csv.pl [-s n,...] [-v] file.xls\n";

$filename =~ /(.*)\.xls$/i or die "Invalid filename: $filename";
my $prefix = $1;

print "Reading '$filename'...\n" if $opt_v;
my $xls = Spreadsheet::ParseExcel::Simple->read($filename)
    or die "Can't read $filename: $!";

my $index = 1;
foreach my $sheet ($xls->sheets) {
    next if $opt_s && !$indexes{$index};

    print "Writing '$prefix-Sheet$index.csv'...\n" if $opt_v;
    open(CSV, ">$prefix-Sheet$index.csv");
    while ($sheet->has_data) {
        my @row = $sheet->next_row;
        foreach (@row) {
            if (/[,"]/) {
                s/"/""/g;
                s/^/"/;
                s/$/"/;
            }
        }
        print CSV join(',', @row), "\n";
    }
    close CSV;
} continue {
    $index++;
}

__END__

=head1 NAME

xls2csv.pl - save MS Excel worksheets as CSV files

=head1 SYNOPSIS

B<xls2csv.pl> [B<-s> n,...] [B<-v>] F<spreadsheet.xls>

=head1 DESCRIPTION

This script saves worksheets from an MS Excel file as CSV files. It does not
depend on Excel because it uses L<Spreadsheet::ParseExcel::Simple> module. The
values in each output CSV file are the underlying data values from each cell in
the Excel file. Thus, cell formatting from the Excel file is not preserved.

The CSV files are saved into the same directory as F<spreadsheet.xls> under the
names F<spreadsheet-SheetI<n>.csv>, where I<n> is the index of each worksheet.

=head1 OPTIONS

=over 7

=item B<-s>

Saves only the worksheets at the specified indexes. More than one index can be
specified with a comma-separated list. The first worksheet has index 1.

=item B<-v>

Be verbose. Some messages indicating the progress will be printed to stdout.

=back

=head1 AUTHOR

Alexander Anderson E<lt>a.anderson@utoronto.caE<gt>

=head1 README

Save MS Excel worksheets as CSV files.

=head1 PREREQUISITES

Spreadsheet::ParseExcel::Simple

=head1 SCRIPT CATEGORIES

Win32
Win32/Utilities

=cut
