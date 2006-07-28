#!/usr/local/bin/perl
# go_runs_tabulate.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

go_runs_tabulate.PLS - DESCRIPTION 

=head1 SYNOPSIS

This is a very shitty script to get tables like this one:

subdir,goid,gostring,orthgenenum,annotgenenum,totalgenenum
hyphymp0001,GO0000118,#histone_deacetylase_complex,1,6,6
hyphymp0002,GO0000119,#mediator_complex,9,24,29
hyphymp0003,GO0000123,#histone_acetyltransferase_complex,2,5,8
hyphymp0004,GO0000124,#SAGA_complex,1,3,4
hyphymp0005,GO0000139,#Golgi_membrane,3,13,14
hyphymp0006,GO0000145,#exocyst,2,8,8
hyphymp0007,GO0000151,#ubiquitin_ligase_complex,20,108,125
hyphymp0008,GO0000152,#nuclear_ubiquitin_ligase_complex,3,14,16
hyphymp0009,GO0000159,#protein_phosphatase_type_2A_complex,1,6,7

Calls will look like:

perl /home/avb/ortholytics/test_lousy_GO_runs_tabulator.PLS -b /home/avb/wallace/eukarya/drosophila/hc_cc_usgo_wDgri_fr_rlong -g /home/avb/wallace/eukarya/drosophila/funccats/nci_nih/gominer_results -c /home/avb/wallace/eukarya/drosophila/_cc_usgo_wDgri_logs -part cc 

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

my ($batchdir,$gcedir,$conclogsdir,
    $outfile,$listfile,$part,$yesnoopt);

GetOptions(
	   'b|batchdir:s' => \$batchdir,
	   'g|gcedir:s' => \$gcedir,
	   'c|conclogsdir:s' => \$conclogsdir,
	   'o|outfile:s' => \$outfile,
	   'l|listfile:s' => \$listfile,
	   'p|part:s' => \$part,
           'yesnoopt' => \$yesnoopt,
          );

my ($subdir, $batchfile);

opendir DIR, $batchdir or die "couldnt find $batchdir:$!";
my ($goid, $gostring,$orthgenenum,$totalgenenum,$annotgenenum,$concatlength,$dummy,@dummyarray);
open OUTFILE, ">$outfile" or die "couldnt open $outfile:$!\n";
open LISTFILE, ">$listfile" or die "couldnt open $listfile:$!\n";
print STDERR "Entries will look like this:\n";
print STDERR "subdir,goid,gostring,orthgenenum,annotgenenum,totalgenenum,concatlength\n";
print OUTFILE "subdir,goid,gostring,orthgenenum,annotgenenum,totalgenenum,concatlength\n";
while (defined($subdir = readdir(DIR))) {
    next if ($subdir =~ /\./);
    next if ($subdir =~ /\.\./);
    opendir SUBDIR, "$batchdir/$subdir" or die "couldnt find $batchdir/$subdir:$!";
    while (( defined($batchfile = readdir(SUBDIR)) )) {
        if ($batchfile =~ /(\S+)\.bf$/ || $batchfile =~ /(\S+)\.ctl$/) {
            open BATCHFILE, "$batchdir/$subdir/$batchfile" or die "couldnt open $batchdir/$subdir/$batchfile";
            my ($goid, $gcefilename);
            while (<BATCHFILE>) {
                if ($_ =~ /.+(GO\d+).fasta/ || /.+(GO\d+).seqt.phy/) {
                    $goid = $1;
                    $gcefilename = "$part.gce.$goid";
                    $gostring = `head -n 2 $gcedir/$gcefilename | tail -n 1`;
                    $gostring =~ s/\n//;
                    $gostring = "\"$gostring\"";
                    ($orthgenenum,$dummy) = split /\ /, `wc -l $conclogsdir/gene_lengths.concatenome.$goid.log`;
                    open GENELENGTHS, "$conclogsdir/gene_lengths.concatenome.$goid.log" 
                        or die "couldnt open $conclogsdir/gene_lengths.concatenome.$goid.log:$!";
                    if ($orthgenenum) {
                        my @genelist;
                        while (<GENELENGTHS>) {
                            if (/.+\;(\S+)\;\S+$/) {
                                push @genelist, $1;
                            } else {
                                print STDERR "parsing error\n";
                            }
                        }
                        close GENELENGTHS;
                        my $listfilestring = "$subdir:";
                        foreach my $gene (@genelist) {
                            $listfilestring .= "$gene,";
                        }
                        $listfilestring =~ s/\,$/\n/;
                        print LISTFILE $listfilestring;
                    }
                    @dummyarray = split /\ /, `grep \'final length\' $conclogsdir/logfile.concatenome.$goid.log`;
                    $concatlength = $dummyarray[-1];
                    $concatlength =~ s/\n//;
                    ($totalgenenum,$dummy) = split /\ /, `wc -l $gcedir/$gcefilename`;
                    $totalgenenum = $totalgenenum - 2;
                    ($annotgenenum,$dummy) = split /\ /, `grep "^#" $gcedir/$gcefilename | wc -l`;
                    $annotgenenum = $totalgenenum - $annotgenenum + 2;
                }
            }
            print OUTFILE "$subdir,$goid,$gostring,$orthgenenum,$annotgenenum,$totalgenenum,$concatlength\n";
        }
    }
}

close OUTFILE;
close LISTFILE;

1;
