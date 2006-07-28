#!/usr/local/bin/perl
# codon_freqs_bias.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

codon_freqs_bias.PLS - Obtain the codon frequencies and codon biases
of the more present amino acids (e.g., >=10 - "n" option) in a gene
from the ENCPrimer codcnt and result files for that genome
CDSs. Calculate the linear regression for each codon freq to the Ncp,
Nc, ScaledChi and B_KM value, and decide the state (P,U or N) for each
codon as a combination of Ncp And ScaledChi states.

This script depends on Statistics::OLS and Statistics::Distributions,
which can be found in CPAN repositories.

States at each codon are:

P - preferred with p-value < 0.01
p - preferred with 0.01 < p-value < 0.05
U - unpreferred with p-value < 0.01
u - preferred with 0.01 < p-value < 0.05
N - neutral (p-value not significant)
n - unique translation - for Trp, Met, etc (varies with codontable)

The NCBI notation corresponds to codons as follows:

Bases at each position are:

 Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
 Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
 Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

=head1 SYNOPSIS

perl codon_freqs_bias.PLS \
-n 20 \
-t 11 \
-c _NC_006624.fa.codcnt \
-e _NC_006624.fa.ENCPrimeResults.main.txt \
[-v]

=head1 DESCRIPTION

-n 20 - minimum number of codons for that aa in that gene to be used in the freqs
-t 11 - translation table (NCBI tables)
-c ENCprime codcnt file for the genome
-e ENCPrime results file for the genome
[-v] verbose - print extra stuff: freq entries, p-values, sample sizes 

=head1 AUTHOR - Albert Vilella

Email 

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut


# Let the code begin...

use strict;
use Getopt::Long;
use Bio::Tools::CodonTable;
use Statistics::OLS;
use Statistics::Distributions;

my ($numseq,
    $tableid,
    $totalcodons,
    $numcod,
    $codcnt_file,
    $encres_file,
    $verbose,
    @codons,
    $codontable,
   ) = "";

my (%genefreqs,
    %genebiases,
    %geneaas,
    %codonfreqs,
    %aaentries,
    %gene_codonentries,
    %pun_table);

GetOptions(
 	   't|table|:s'=> \$tableid,
 	   'n|num|numcod:s'=> \$numcod,
	   'c|codcnt:s' => \$codcnt_file,
	   'e|enc|encres:s' => \$encres_file,
           'v|verbose' => \$verbose,
          );

$codontable = Bio::Tools::CodonTable->new( -id => $tableid);

open (CODCNT,"$codcnt_file") or die "cannot open file $codcnt_file: $!";
print "Loading codon counts...\n";
# 1897
# 64
$numseq = <CODCNT>;
$totalcodons = <CODCNT>;
# TTT TTC TTA TTG TCT TCC TCA TCG TAT TAC TAA TAG TGT TGC TGA TGG CTT CTC CTA CTG CCT CCC CCA CCG CAT CAC CAA CAG CGT CGC CGA CGG ATT ATC ATA ATG ACT ACC ACA ACG AAT AAC AAA AAG AGT AGC AGA AGG GTT GTC GTA GTG GCT GCC GCA GCG GAT GAC GAA GAG GGT GGC GGA GGG 
my $line = <CODCNT>;
@codons = ($line =~ /\S+/gi);
while (<CODCNT>) {
    # 1452198>0 2 5 8 1 2 5 10 0 1 6 7 2 4 7 7 0 3 5 2 1 3 4 5 1 1 2 7 0 1 4 6 0 0 2 4 4 4 4 5 0 3 5 6 5 4 19 9 0 3 6 3 0 2 1 1 0 1 3 2 3 3 4 2 
    #...
    #Totals> 9065 14870 9234 8153 3803 3947 4944 2886 7751 13264 559 396 1367 1679 960 6509 11742 11389 10177 5844 5040 5446 9228 3564 3138 5080 4010 5103 570 660 497 384 11738 8478 26460 12759 6096 5979 4455 6518 7165 11146 13541 29293 4112 7664 9307 20087 21448 8210 7476 7506 11044 10562 9209 5654 15110 10134 19597 28953 8451 5990 16474 8814 
     next if ($_ =~ /^Totals/);
     if ($_ =~ /(\S+)\>(.+)/g) {
         my $gene = $1;
         my @freqs = ($2 =~ /\S+/gi);
         my $i = 0;
         for ($i=0; $i<$totalcodons; $i++) {
             $genefreqs{$gene}{$codons[$i]} = $freqs[$i];
         }
     } else {
         warn "parsing error: $!";
     }
}
close(CODCNT);

open (ENCRES,"$encres_file") or die "cannot open file $encres_file: $!";
print "Loading ENC Results...\n";
while (<ENCRES>) {
     #Name;Nc;Ncp;Sch;SumChi;df;p;BKM;n_codons
     #1452198:;52,4766;54,6991;0,2650;58,2913;43;0,0598;0,3771;220
     #...
     #Totals:;49,7731;50,8532;0,1843;101515,1839;43;0,0000;0,3417;550689
     next if ($_ =~ /^Totals/);
     next if ($_ =~ /^Name/);
     if ($_ =~ /^(\S+)\:\;(\S+)\;(\S+)\;(\S+)\;(\S+)\;(\S+)\;(\S+)\;(\S+)\;(\S+)/
         ||
         $_ =~ /^(\S+)\:\s(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)/) {
         my $gene = $1;
         $genebiases{$gene}{Nc} =          $2;
         $genebiases{$gene}{Ncp} =         $3;
         $genebiases{$gene}{Sch} =         $4;
         $genebiases{$gene}{SumChi} =      $5;
         $genebiases{$gene}{df} =          $6;
         $genebiases{$gene}{p} =           $7;
         $genebiases{$gene}{BKM} =         $8;
         $genebiases{$gene}{n_codons} =    $9;
     } else {
         warn "parsing error: $!";
     }
}
close(ENCRES);

print "Counting codons...\n";
foreach my $gene (keys %genefreqs) {
    foreach my $codon (keys %{ $genefreqs{$gene} }) {
        my $aa = $codontable->translate($codon);
        $geneaas{$gene}{$aa}{$codon} = $genefreqs{$gene}{$codon};
        $geneaas{$gene}{$aa}{'_codcnt'} += $genefreqs{$gene}{$codon};
    }
}

$numcod = sprintf("%02d",$numcod);
print "Parsing chosen entries...\n";
foreach my $gene (keys %geneaas) {
    foreach my $aa (keys %{ $geneaas{$gene} }) {
        if ($geneaas{$gene}{$aa}{'_codcnt'} >= $numcod) {
            foreach my $codon (sort keys %{ $geneaas{$gene}{$aa}} ) {
                my $freq = $geneaas{$gene}{$aa}{$codon}/$geneaas{$gene}{$aa}{'_codcnt'};
                $freq = sprintf("%05f", $freq);
                unless ($codon =~ /codcnt/) {
                    my $aastring = ""; my $entry = "";
                    $entry = "$aa" . "_" . "$codon" . "_". "$gene";
                    $aastring = "$aa $codon $gene $geneaas{$gene}{$aa}{$codon} $freq";
                    $aastring .= " ";
                    $aastring .= "$genebiases{$gene}{Ncp} $genebiases{$gene}{Nc} $genebiases{$gene}{Sch} $genebiases{$gene}{SumChi} $genebiases{$gene}{df} $genebiases{$gene}{p} $genebiases{$gene}{BKM} $genebiases{$gene}{n_codons}";
                    $aaentries{$entry} = $aastring;
                    # codon entries for linear regressions
                    $gene_codonentries{$codon}{Ncp} .= "$genebiases{$gene}{Ncp} $freq ";
                    $gene_codonentries{$codon}{Nc} .= "$genebiases{$gene}{Nc} $freq ";
                    $gene_codonentries{$codon}{Sch} .= "$genebiases{$gene}{Sch} $freq ";
                    $gene_codonentries{$codon}{BKM} .= "$genebiases{$gene}{BKM} $freq ";
                }
            }
        }
    }
}

if ($verbose) {
    print "Writing freqs.bias outfile $encres_file.$numcod.freqs.bias.txt ...\n";
    open(OUTFILE, ">$encres_file.$numcod.freqs.bias.txt") or die $!;
    print OUTFILE "gene aa codon codcnt codfreq Ncp Nc Sch SumChi df p BKM n_codons\n";
    foreach my $entry (sort keys %aaentries) {
        print OUTFILE "$aaentries{$entry}\n";
    }
    close(OUTFILE);
}

# ENC
# Ncp
# Sch
# BKM

print "Calculating regressions...\n";
print "Writing REGRFILE $encres_file.$numcod.regr.txt ...\n";
open(REGRFILE, ">$encres_file.$numcod.regr.txt") or die $!;
my ($Ncp_signif_pref, $ENC_signif_pref, $Sch_signif_pref, $BKM_signif_pref,
    $Ncp_signif_unpref, $ENC_signif_unpref, $Sch_signif_unpref, $BKM_signif_unpref
   ) = 0;
my (@Ncp_preferred, @ENC_preferred, @Sch_preferred, @BKM_preferred) = "";
foreach my $codon (sort keys %gene_codonentries) {
    my $aa = $codontable->translate($codon);
    my @Ncp_xydataset = split(/\ /, $gene_codonentries{$codon}{Ncp});
    my @ENC_xydataset = split(/\ /, $gene_codonentries{$codon}{Nc});
    my @Sch_xydataset = split(/\ /, $gene_codonentries{$codon}{Sch});
    my @BKM_xydataset = split(/\ /, $gene_codonentries{$codon}{BKM});

    if (2 >= scalar(@Ncp_xydataset)) {
        $pun_table{$aa}{$codon}{'sample_size'} = 0;
    } else {
        my $ls_Ncp = Statistics::OLS->new();
        my $ls_ENC = Statistics::OLS->new();
        my $ls_Sch = Statistics::OLS->new();
        my $ls_BKM = Statistics::OLS->new();

        $ls_Ncp->setData (\@Ncp_xydataset);
        $ls_ENC->setData (\@ENC_xydataset);
        $ls_Sch->setData (\@Sch_xydataset);
        $ls_BKM->setData (\@BKM_xydataset);

        $ls_Ncp->regress();
        $ls_ENC->regress();
        $ls_Sch->regress();
        $ls_BKM->regress();

        my $sample_size = $ls_Ncp->size();
        $sample_size = sprintf("%05d",$sample_size);
        $pun_table{$aa}{'sample_size'} = $sample_size;

        my ($Ncp_intercept, $Ncp_slope) = $ls_Ncp->coefficients();
        my ($ENC_intercept, $ENC_slope) = $ls_ENC->coefficients();
        my ($Sch_intercept, $Sch_slope) = $ls_Sch->coefficients();
        my ($BKM_intercept, $BKM_slope) = $ls_BKM->coefficients();

        my ($Ncp_tstat_intercept, $Ncp_tstat_slope) = $ls_Ncp->tstats();
        my ($ENC_tstat_intercept, $ENC_tstat_slope) = $ls_ENC->tstats();
        my ($Sch_tstat_intercept, $Sch_tstat_slope) = $ls_Sch->tstats();
        my ($BKM_tstat_intercept, $BKM_tstat_slope) = $ls_BKM->tstats();

        my $Ncp_tprob = 2*Statistics::Distributions::tprob((($sample_size)-2),abs($Ncp_tstat_slope));
        my $ENC_tprob = 2*Statistics::Distributions::tprob((($sample_size)-2),abs($ENC_tstat_slope));
        my $Sch_tprob = 2*Statistics::Distributions::tprob((($sample_size)-2),abs($Sch_tstat_slope));
        my $BKM_tprob = 2*Statistics::Distributions::tprob((($sample_size)-2),abs($BKM_tstat_slope));

        $pun_table{$aa}{$codon}{'Ncp_method'}{'pvalue'} = $Ncp_tprob;
        $pun_table{$aa}{$codon}{'ENC_method'}{'pvalue'} = $ENC_tprob;
        $pun_table{$aa}{$codon}{'Sch_method'}{'pvalue'} = $Sch_tprob;
        $pun_table{$aa}{$codon}{'BKM_method'}{'pvalue'} = $BKM_tprob;

        $pun_table{$aa}{$codon}{'Ncp_method'}{'slope'} = $Ncp_tstat_slope;
        $pun_table{$aa}{$codon}{'ENC_method'}{'slope'} = $ENC_tstat_slope;
        #Slopes for these two methods are modified so that all have to
        #go to the same direction
        $pun_table{$aa}{$codon}{'Sch_method'}{'slope'} = -($Sch_tstat_slope);
        $pun_table{$aa}{$codon}{'BKM_method'}{'slope'} = -($BKM_tstat_slope);
    }
}

# Pref/Unpref significance
foreach my $aa (keys %pun_table) {
    foreach my $codon (keys %{ $pun_table{$aa} }) {
        next if $codon =~ /sample_size/;
        foreach my $method (keys %{ $pun_table{$aa}{$codon} }) {
            next unless $method =~ /method/;
            if (    ($pun_table{$aa}{$codon}{$method}{'pvalue'} < 0.01) && ($pun_table{$aa}{$codon}{$method}{'slope'} < 1) )
                {$pun_table{$aa}{$codon}{$method}{'pun'} = 'P';}
            elsif ( ($pun_table{$aa}{$codon}{$method}{'pvalue'} < 0.05) && ($pun_table{$aa}{$codon}{$method}{'slope'} < 1) )
                {$pun_table{$aa}{$codon}{$method}{'pun'} = 'p';}
            elsif ( ($pun_table{$aa}{$codon}{$method}{'pvalue'} < 0.01) && ($pun_table{$aa}{$codon}{$method}{'slope'} > 1) )
                {$pun_table{$aa}{$codon}{$method}{'pun'} = 'U';}
            elsif ( ($pun_table{$aa}{$codon}{$method}{'pvalue'} < 0.05) && ($pun_table{$aa}{$codon}{$method}{'slope'} > 1) )
                {$pun_table{$aa}{$codon}{$method}{'pun'} = 'u';}
            else 
                {$pun_table{$aa}{$codon}{$method}{'pun'} = 'N';}
        }
    }
}

# Decision making
# Ncp + Sch => ok
# u    + U   => u
# p    + P   => p
# p    + u   => N
# p,P  + N   => N
# u,U  + N   => N
foreach my $aa (sort keys %pun_table) {
    next if $aa =~ /\*/;
    foreach my $codon (sort keys %{ $pun_table{$aa} }) {
        next if $codon =~ /sample_size/;
        my $a = $pun_table{$aa}{$codon}{'Ncp_method'}{'pun'};
        $a .= $pun_table{$aa}{$codon}{'Sch_method'}{'pun'};
        if (
            ($a =~ /N/) ) {
            $pun_table{$aa}{$codon}{'result'}{'pun'} = 'N';
        } elsif (
            ($a =~ /UU/) ) {
            $pun_table{$aa}{$codon}{'result'}{'pun'} = 'U';
        } elsif (
            ($a =~ /PP/) ) {
            $pun_table{$aa}{$codon}{'result'}{'pun'} = 'P';
        } elsif (
            ($a =~ /uu/) ) {
            $pun_table{$aa}{$codon}{'result'}{'pun'} = 'u';
        } elsif (
            ($a =~ /pp/) ) {
            $pun_table{$aa}{$codon}{'result'}{'pun'} = 'p';
        } elsif (
            ($a =~ /Pp/) ) {
            $pun_table{$aa}{$codon}{'result'}{'pun'} = 'p';
        } elsif (
            ($a =~ /pP/) ) {
            $pun_table{$aa}{$codon}{'result'}{'pun'} = 'p';
        } elsif (
            ($a =~ /uU/) ) {
            $pun_table{$aa}{$codon}{'result'}{'pun'} = 'u';
        } elsif (
            ($a =~ /Uu/) ) {
            $pun_table{$aa}{$codon}{'result'}{'pun'} = 'u';
        } elsif (
            ($a =~ /UP/) ) {
            $pun_table{$aa}{$codon}{'result'}{'pun'} = 'N';
        } elsif (
            ($a =~ /PU/) ) {
            $pun_table{$aa}{$codon}{'result'}{'pun'} = 'N';
        } elsif (
            ($a =~ /up/) ) {
            $pun_table{$aa}{$codon}{'result'}{'pun'} = 'N';
        } elsif (
            ($a =~ /pu/) ) {
            $pun_table{$aa}{$codon}{'result'}{'pun'} = 'N';
        } elsif (
            ($a =~ /Up/) ) {
            $pun_table{$aa}{$codon}{'result'}{'pun'} = 'N';
        } elsif (
            ($a =~ /pU/) ) {
            $pun_table{$aa}{$codon}{'result'}{'pun'} = 'N';
        } elsif (
            ($a =~ /uP/) ) {
            $pun_table{$aa}{$codon}{'result'}{'pun'} = 'N';
        } elsif (
            ($a =~ /Pu/) ) {
            $pun_table{$aa}{$codon}{'result'}{'pun'} = 'N';
        } else {
            warn "error in decision making\n";
        }
    }
}

#Hash in NCBI order
my %pun_table_codons;
my @nucl_array = qw[T C A G];
foreach my $base1 (@nucl_array) {
    foreach my $base2 (@nucl_array) {
        foreach my $base3 (@nucl_array) {
            my $codon = $base1 . $base2 . $base3;
            $pun_table_codons{$codon} = 'N';
        }
    }
}
foreach my $aa (sort keys %pun_table) {
    next if $aa =~ /\*/;
    foreach my $codon (sort keys %{ $pun_table{$aa} }) {
        next if $codon =~ /sample_size/;
        $pun_table_codons{$codon} = $pun_table{$aa}{$codon}{'result'}{'pun'};
        #Hard code non-variable AAs as 'n'
        if (1 == scalar($codontable->revtranslate($aa))) {
            $pun_table_codons{$codon} = 'n';
            $pun_table{$aa}{$codon}{'result'}{'pun'} = 'n';
        }
    }
}

# Printing
foreach my $aa (sort keys %pun_table) {
    my $string;
    $string .= "("."$aa".")"."$pun_table{$aa}{sample_size}"."\n";
    print REGRFILE $string;print $string;$string = "";
    foreach my $codon (sort keys %{ $pun_table{$aa} }) {
        next if $codon =~ /sample_size/;
        $string .= "("."$aa".")"."$codon ";
        $string .= "$pun_table{$aa}{$codon}{'result'}{'pun'}";
        print REGRFILE "$string\n";print "$string\n";$string = "";
    }
}

if ($verbose) {
    foreach my $aa (sort keys %pun_table) {
        my $string;
        $string .= "("."$aa".")"."$pun_table{$aa}{sample_size}"."\n";
        print REGRFILE $string;print $string;$string = "";
        foreach my $codon (sort keys %{ $pun_table{$aa} }) {
            next if $codon =~ /sample_size/;
            $string .= "("."$aa".")"."$codon ";
            foreach my $method (sort keys %{ $pun_table{$aa}{$codon} }) {
                next unless $method =~ /method/;
                $string .= "$pun_table{$aa}{$codon}{$method}{'pvalue'}"."($pun_table{$aa}{$codon}{$method}{'slope'})"." ";
            }
            print REGRFILE "$string\n";print "$string\n";$string = "";
        }
    }
}

print REGRFILE "\n";print "\n";

#Printing in NCBI order
my $string = "";
my $count = 1;
foreach my $base1 (@nucl_array) {
    foreach my $base2 (@nucl_array) {
        foreach my $base3 (@nucl_array) {
            my $codon = $base1 . $base2 . $base3;
            $string .= "$pun_table_codons{$codon}";
        }
    }
}
print REGRFILE "$string\n";print "$string\n";$string = "";

close(REGRFILE);

1;
