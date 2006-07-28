#!/usr/local/bin/perl

# cog_histograms.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

cog_histograms.PLS - DESCRIPTION 

=head1 SYNOPSIS

perl cog_histograms.PLS NC_0000001.ptt

=head1 DESCRIPTION

Count the number of cogs and give an histogram of the number of genes
in each cog identifier.

Options:

ptt files - ptt extension

See NCBI ptt documentation

trio files - trio extension

Showing only tuples.
Field separator is "	".
 601.GSA_SALTI.AL513382.235499   | 602.GSA_SALTY.AE006468.238491   | 83333.GSA_ECOLI.U00096.174882   | COG0001             | COG0001             | COG0001               |                          235499 |                          238491 |                            174882 |                        234219 |                        237211 |                          173602 |               601 |               602 |               83333 |                 426 |                 426 |                   426


=head1 AUTHOR - Albert Vilella

Email 

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut


# Let the code begin...

use strict;

1;

use strict;
use POSIX;          #Used for floor
use Carp;           #Used for read_file and write_file
use Fcntl;          #Used for read_file and write_file
use File::Basename; #Used for filename and dirname parsing
use File::Spec;     #Used for filename and dirname parsing
use Getopt::Long;   #Used for GetOptions

use Bio::Graphics;
#use GD::Graph::histogram;
#use Bio::Graphics::Glyph::xyplot;


my $inputfile = shift;
my %inputfile;
($inputfile{'base'}, $inputfile{'dir'}, $inputfile{'ext'}) = fileparse($inputfile,'\.\w+?');

my $cog;
my $bin_num = 0;
my $total;
my $withcog;
my $nocog;
my $description;
my %cog_hash;
my %occurrences_hash;
my %histogram_entries;

open INPUTFILE,"$inputfile" or die "$!";

print "Reading input file\n";
if ($inputfile{'ext'} =~ /trio/) {
    while (<INPUTFILE>) {
        next if ($_ =~ /^Showing only tuples./);
        next if ($_ =~ /^Field separator is "	"./);
        next if ($_ =~ /^\s*$/);         # (
        if ($_ =~ /\s*\S+\s+\|\s*\S+\s+\|\s*\S+\s+\|\s*(\SOG\S+)\s+\|/) { 
            if (($1 == $2) && ($2 == $3)) {
                $cog_hash{$1}++;
                $withcog++;
                $total++;
            } else {
                print "Mismatched cogs: $1\t$2\t$3 (withcog $withcog ; total $total)\n";
            }
        } else {
            die "Line mismatch: (withcog $withcog ; total $total) $!\n";
        }
    }
} elsif ($inputfile{'ext'} =~ /ptt/) {
    while (<INPUTFILE>) {
        next if ($_ =~ /^Location/);
        if (     $_ =~ /\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t(COG\S+)\s+/) {
            $cog_hash{$1}++;
            $withcog++;
            $total++;
        } elsif ($_ =~ /\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t(\-*)\s+/) {
            $nocog++;
            $total++;
        }
    }
} else {
    print "Extension of input file is unknown. Read documentation.\n";
    exit(1);
}

close INPUTFILE;

if ($withcog == 0) {
    print "Sorry. No COGs in this file!\n";
    exit;
}

print "Binning (this might take a while)\n";
$histogram_entries{0} = 0;
my $total_occurrences = 0;
foreach $cog (keys %cog_hash) {
    $total_occurrences = $cog_hash{$cog};
    $total_occurrences = sprintf("%07d", $total_occurrences);
    $occurrences_hash{$total_occurrences} .= "$cog;";
    $histogram_entries{$total_occurrences}++;
}

$bin_num = scalar(keys %histogram_entries);
my $total_cogs = scalar(keys %cog_hash);

my %outfile;
$outfile{'dir'} = File::Spec->curdir();
$outfile{'base'} = $inputfile{'base'};
$outfile{'ext'} = 'bed';

my $outfile = File::Spec->catfile( $outfile{'dir'}, "$outfile{'base'}\.$outfile{'ext'}" );

my $triangular_cogs;
open(OUTFILE, ">$inputfile.percent.txt") or die $!;
foreach my $entry (sort keys %histogram_entries) {
    next if ($entry == 0);
    if ($entry == 1) {
        $triangular_cogs = $histogram_entries{$entry};
        my $percent = sprintf("%2.2f",($triangular_cogs*100 /$total_cogs));
        $description = "$percent% ($triangular_cogs/$total_cogs) of triangular COGs out of $total_cogs COGs for $inputfile{'base'}";
        print OUTFILE "$description\n";
    }
    print OUTFILE "$entry\t$histogram_entries{$entry}\n";
}

print OUTFILE "\nCOGs combinations occurrence\n";
foreach $cog (sort keys %cog_hash) {
    print OUTFILE "$cog\t$cog_hash{$cog}\n";
}

print OUTFILE "\noccurrences COG combinations\n";
foreach my $occurrence (sort keys %occurrences_hash) {
    print OUTFILE "$occurrence\t$occurrences_hash{$occurrence}\n\n";
}


print "Printing text file\n";
close OUTFILE;


my $ftr = 'Bio::Graphics::Feature';
my $segment;
my $width = 820;
my $length = 1000;
my $pad_left = 10;
my $pad_right = 10;
my $start = 0;
my $fontsize = 170;
my $segment_width = floor($width/$bin_num);
my $end = $segment_width;

my $panel = Bio::Graphics::Panel->new(
                                      -length => "$length",
                                      -height => "$length",
                                      -width  => "$width",
                                      -pad_left => "$pad_left",
                                      -pad_right => "$pad_right",
                                     );

my $svgpanel = Bio::Graphics::Panel->new(
                                         -image_class=>'SVG',
                                         -length => "$length",
                                         -height => "$length",
                                         -width  => "$width",
                                         -pad_left => "$pad_left",
                                         -pad_right => "$pad_right",
                                        );


my $partial_gene = $ftr->new(-name   =>"$inputfile{'base'}",
			     -strand => '+1',
			     -type   =>"$inputfile{'ext'} file",
			     -source =>"$inputfile{'base'}");

my $score=0;
foreach my $entry (sort keys %histogram_entries) {
    next if ($entry == 0);
    if ($entry == 1) {
        $triangular_cogs = $histogram_entries{$entry};
        my $percent = sprintf("%2.2f",($triangular_cogs*100 /$total_cogs));
        $description = "$percent% ($triangular_cogs/$total_cogs) of triangular COGs out of $total_cogs COGs for $inputfile{'base'}";
    }
    $start = $start + $segment_width;
    $end = $end + $segment_width;
    $score = $histogram_entries{$entry};
    $segment = $ftr->new(-start=>$start,
                         -end=>$end,
                         -key=> "$entry",
                         -name=>" $entry",
                         -score=>"$score"
                        );
    $partial_gene->add_segment($segment);
}

my $track = $panel->add_track($partial_gene,
                               -glyph => 'xyplot',
                              -length => "$length",
                              -height => "$length",
                              -label  => 1,
                              -bgcolor => 'white',
                              -fgcolor => 'blue',
                              -font2color => 'black',
                              -fontsize => "$fontsize",
                              -description => "$description",
                              -graph_type => 'boxes',
                             );

my $svgtrack = $svgpanel->add_track($partial_gene,
                                    -glyph => 'xyplot',
                                    -length => "$length",
                                    -height => "$length",
                                    -label  => 1,
                                    -bgcolor => 'white',
                                    -fgcolor => 'blue',
                                    -font2color => 'black',
                                    -fontsize => "$fontsize",
                                    -description => "$description",
                                    -graph_type => 'boxes',
                                   );


#Make each segment to be like a partial gene
open(OUTFILE, ">$inputfile.percent.png") or die $!;
binmode OUTFILE;
print "Printing graphic file\n";
print OUTFILE  $panel->png;
close OUTFILE;

open(OUTFILE, ">$inputfile.percent.svg") or die $!;
binmode OUTFILE;
print "Printing graphic file\n";
print OUTFILE  $svgpanel->svg;
close OUTFILE;


#my $graph_percent = new GD::Graph::histogram(800,600);

# $graph_percent->set( 
#             x_label         => 'number of genes per COG group', #'
#             y_label         => 'Percentage',
#             title           => "$description";
#             x_labels_vertical => 1,
#             bar_spacing     => 0,
#             shadow_depth    => 1,
#             shadowclr       => 'dred',
#             transparent     => 0,
#             histogram_type  => 'percentage',
#             histogram_bins  => $bin_num,
#            ) 
#     or warn $graph_percent->error;

# print "Flushing data to plot object\n";
#my $gd_percent = $graph_percent->plot(\@data) or die $graph_percent->error;

# open(OUTFILE, ">$inputfile.percent.png") or die $!;
# binmode OUTFILE;
# print "Printing graphic file\n";
# print OUTFILE $gd_percent->png;
# close OUTFILE;

# my $graph_count = new GD::Graph::histogram(800,600);

# $graph_count->set( 
#             x_label         => 'X Label', #'
#             y_label         => 'Count',
#             title           => "$description",
#             x_labels_vertical => 1,
#             bar_spacing     => 0,
#             shadow_depth    => 1,
#             shadowclr       => 'dred',
#             transparent     => 0,
#             histogram_type  => 'count',
#             histogram_bins  => $bin_num,
#            ) 
#     or warn $graph_count->error;

# print "Flushing data to plot object\n";
# my $gd_count = $graph_count->plot(\@data) or die $graph_count->error;

# open(OUTFILE, ">$inputfile.count.png") or die $!;
# print "Output file $outfile\n";
# binmode OUTFILE;
# print OUTFILE $gd_count->png;
# print "Printing graphic file\n";
# close OUTFILE;


1;
