# calculate_prefunpref.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

calculate_prefunpref.PLS - Calculate the changes of state between
preferred and unpreferred codons between a pair of aligned coding
sequences

=head1 SYNOPSIS

perl calculate_prefunpref.PLS
-aln myfile.fasta
-table UPUNUNUNUPNNNNNnNPNNUpNNUPUPUUUUUNPnUPUPUPUPUPUPPPUUNPUNUPUPNNNu
[-table2 PPNNNNUNUPNNNNNnNPNNUpNNUPUPUUUUUNPnUPUPUPUPUPUPPPUUNPUNUPUPNNNu]
[-simplified]

=head1 DESCRIPTION

[-table2 PPNNNNUNUPNNNNNnNPNNUpNNUPUPUUUUUNPnUPUPUPUPUPUPPPUUNPUNUPUPNNNu]
    A second table is used with this option, then a common set from
    both tables is created

[-simplified]
    Calculate also the simplified pref/nopref [P/A] from Akashi.

[-uppercase]
    Consider all "u" and "p" as uppercased-values (signif 0.01 or 0.05).

[-outxls /my/path/myfile.xls]
    Generate an xls file with the tables (requires Spreadsheet::WriteExcel)

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


=head1 AUTHOR - Albert Vilella

Email 

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut

# Let the code begin...

use strict;
use Getopt::Long;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Tools::GuessSeqFormat;
use Bio::Align::DNAStatistics;

my ($alnfile, $puntable1, $puntable2,
    $simplified, $uppercase, $outxls) = undef;

GetOptions(
 	   'a|aln|align|alignment:s'=> \$alnfile,
	   't|t1|tab|table|table1:s' => \$puntable1,
	   't2|table2:s' => \$puntable2,
           's|sim|simpl|simplified' => \$simplified,
           'u|uc|upper|uppercase' => \$uppercase,
           'o|outxls:s' => \$outxls,
          );

if ($uppercase) {
    print "Considering all uppercase...\n";
    $puntable1 = uc($puntable1);
    print "$puntable1\n";
    if ($puntable2) {
        $puntable2 = uc($puntable2);
        print "$puntable2\n";
    }
}

my $puntable_strict = undef;
if ($puntable2) {
    $puntable_strict = "";
    my ($equal,$different)= 0;
    print "Defining common sets from 2 tables...\n";
    die "tables of different length$!\n" if ((length $puntable1) != (length $puntable2));
    foreach my $index (0..(length($puntable1)-1)) {
        my $char1 = substr($puntable1, $index, 1);
        my $char2 = substr($puntable2, $index, 1);
        my $char;
        if (($char1) eq ($char2)) {
            $char = $char1;
            $equal++;
        } elsif ((!($char1 eq $char2)) && (uc($char1) eq uc($char2)) ) {
            $char = uc($char1);
            $equal++;
        } else {
            $char = 'x';
            $different++;
        }
        $puntable_strict .= $char;
    }
    print "Equal = $equal; Different = $different\n";
    print "Marked different states with \"x\"\n";
    print "$puntable_strict\n";
}

my $guessed_format = new Bio::Tools::GuessSeqFormat
    (-file=>"$alnfile")->guess;

my $alignin = Bio::AlignIO->new
    (
     -file=>"$alnfile",
     -format=>$guessed_format,
    );

print "Reading pairwise alignment...\n";
my $alnobj = $alignin->next_aln;
my ($seq1id,$seq2id) = map { $_->display_id } $alnobj->each_seq;
print "seq1 is $seq1id\n";
print "seq2 is $seq2id\n";

my @seqs;
for my $seq ($alnobj->each_seq) {
    push @seqs, {id => $seq->display_id, seq=>$seq->seq};
#     ########################################
#     my $debugseq = $seq->subseq(1,90000);
#     push @seqs, {id => $seq->display_id, seq=>$debugseq};
#     ########################################
}

print "Starting DNAStatistics...\n";
my $stats = new Bio::Align::DNAStatistics;

my $seq1 = $seqs[0]{'seq'};
my $seq2 = $seqs[1]{'seq'};

my %synchanges = $stats->get_syn_changes();
my $CODONS = $stats->get_codons();

my %simpl_results = undef;
if ($simplified) {
    print "Simplified prefunpref...\n";
    my $simpl_puntable = $puntable1;
    $simpl_puntable =~ tr/[UuNn]/A/;
    %simpl_results = count_prefunpref
        (
         $seq1, $seq2, $simpl_puntable, $CODONS, \%synchanges
        );

    print "pref/unpref from $seq1id to $seq2id\n";
    for my $key (sort keys %simpl_results) {
        print "$key = $simpl_results{$key}\n";
    }
}

print "Analysis prefunpref...\n";
my %results = count_prefunpref
    (
     $seq1, $seq2, $puntable1, $CODONS, \%synchanges
    );

print "pref/unpref from $seq1id to $seq2id\n";
for my $key (sort keys %results) {
    print "$key = $results{$key}\n";
}

my %results_strict;
my %simpl_results_strict;
if ($puntable2) {
    if ($simplified) {
        print "Strict Simplified prefunpref...\n";
        my $simpl_puntable_strict = $puntable_strict;
        $simpl_puntable_strict =~ tr/[UuNn]/A/;
        %simpl_results_strict = count_prefunpref
            (
             $seq1, $seq2, $simpl_puntable_strict, $CODONS, \%synchanges
            );

        print "pref/unpref from $seq1id to $seq2id\n";
        for my $key (sort keys %simpl_results_strict) {
            print "$key = $simpl_results_strict{$key}\n";
        }
    }

    print "Strict Analysis prefunpref...\n";
    %results_strict = count_prefunpref
        (
         $seq1, $seq2, $puntable_strict, $CODONS, \%synchanges
        );

    print "pref/unpref from $seq1id to $seq2id\n";
    for my $key (sort keys %results_strict) {
        print "$key = $results_strict{$key}\n";
    }
}

################################################################################

# Excel tables
if ($outxls) {
    print "Generating xls file...\n";
    eval { require Spreadsheet::WriteExcel; };
	if ( $@ ) {
            print STDERR "Spreadsheet::WriteExcel not found in the system. Wont generate xls file\n";
        }

    my $workbook = Spreadsheet::WriteExcel->new("$outxls");
    my $worksheet = $workbook->add_worksheet("PUN tables");
    $worksheet->activate();
    $worksheet->set_column(5, 5, 2);

    my $boldrightformat = $workbook->add_format();
    $boldrightformat->set_bold();
    $boldrightformat->set_align('right');
    $worksheet->write('B1', "$seq1id->$seq2id", $boldrightformat);

    my @header = ("Pref", "Neu", "Unp");
    $worksheet->write_row('C2', \@header, $boldrightformat);
    $worksheet->write_col('B3', \@header, $boldrightformat);

    my @pun_table = 
        (
         [$results{'11_PP'}, $results{'11_PN'}, $results{'11_PU'}],
         [$results{'11_NP'}, $results{'11_NN'}, $results{'11_NU'}],
         [$results{'11_UP'}, $results{'11_UN'}, $results{'11_UU'}]
        );
    $worksheet->write_col('C3', \@pun_table);

    if ($simplified) {
        my $P_noP = $results{'11_PN'} + $results{'11_PU'};
        my $noP_P = $results{'11_NP'} + $results{'11_UP'};
        my $U_noU = $results{'11_UN'} + $results{'11_UP'};
        my $noU_U = $results{'11_NU'} + $results{'11_PU'};
        $worksheet->write('G2', "P->noP", $boldrightformat);
        $worksheet->write('G3', "noP->P", $boldrightformat);
        $worksheet->write('H2', $P_noP);
        $worksheet->write('H3', $noP_P);
        my $P_ChiSquare = sprintf("%.3f", ($P_noP-(($P_noP+$noP_P)/2))*($P_noP-(($P_noP+$noP_P)/2))/(($P_noP+$noP_P)/2)+($noP_P-(($P_noP+$noP_P)/2))*($noP_P-(($P_noP+$noP_P)/2))/(($P_noP+$noP_P)/2) );
        $worksheet->write('G4', "ChiSq", $boldrightformat);
        $worksheet->write('H4', $P_ChiSquare);
        $worksheet->write('I2', "U->noU", $boldrightformat);
        $worksheet->write('I3', "noU->U", $boldrightformat);
        $worksheet->write('J2', $U_noU);
        $worksheet->write('J3', $noU_U);
        my $U_ChiSquare = sprintf("%.3f", ($U_noU-(($U_noU+$noU_U)/2))*($U_noU-(($U_noU+$noU_U)/2))/(($U_noU+$noU_U)/2)+($noU_U-(($U_noU+$noU_U)/2))*($noU_U-(($U_noU+$noU_U)/2))/(($U_noU+$noU_U)/2) );
        $worksheet->write('I4', "ChiSq", $boldrightformat);
        $worksheet->write('J4', $U_ChiSquare);
    }

    if ($puntable2) {
        my @header = ("Pref", "Neu", "Unp");
        $worksheet->write('B8', '(strict)', $boldrightformat);
        $worksheet->write_row('C8', \@header, $boldrightformat);
        $worksheet->write_col('B9', \@header, $boldrightformat);

        my @pun_table = 
            (
             [$results_strict{'11_PP'}, $results_strict{'11_PN'}, $results_strict{'11_PU'}],
             [$results_strict{'11_NP'}, $results_strict{'11_NN'}, $results_strict{'11_NU'}],
             [$results_strict{'11_UP'}, $results_strict{'11_UN'}, $results_strict{'11_UU'}]
            );
        $worksheet->write_col('C9', \@pun_table);

        if ($simplified) {
            my $P_noP = $results_strict{'11_PN'} + $results_strict{'11_PU'};
            my $noP_P = $results_strict{'11_NP'} + $results_strict{'11_UP'};
            my $U_noU = $results_strict{'11_UN'} + $results_strict{'11_UP'};
            my $noU_U = $results_strict{'11_NU'} + $results_strict{'11_PU'};
            $worksheet->write('G8', "P->noP", $boldrightformat);
            $worksheet->write('G9', "noP->P", $boldrightformat);
            $worksheet->write('H8', $P_noP);
            $worksheet->write('H9', $noP_P);
            my $P_ChiSquare = sprintf("%.3f", ($P_noP-(($P_noP+$noP_P)/2))*($P_noP-(($P_noP+$noP_P)/2))/(($P_noP+$noP_P)/2)+($noP_P-(($P_noP+$noP_P)/2))*($noP_P-(($P_noP+$noP_P)/2))/(($P_noP+$noP_P)/2) );
            $worksheet->write('G10', "ChiSq", $boldrightformat);
            $worksheet->write('H10', $P_ChiSquare);
            $worksheet->write('I8', "U->noU", $boldrightformat);
            $worksheet->write('I9', "noU->U", $boldrightformat);
            $worksheet->write('J8', $U_noU);
            $worksheet->write('J9', $noU_U);
            my $U_ChiSquare = sprintf("%.3f", ($U_noU-(($U_noU+$noU_U)/2))*($U_noU-(($U_noU+$noU_U)/2))/(($U_noU+$noU_U)/2)+($noU_U-(($U_noU+$noU_U)/2))*($noU_U-(($U_noU+$noU_U)/2))/(($U_noU+$noU_U)/2) );
            $worksheet->write('I10', "ChiSq", $boldrightformat);
            $worksheet->write('J10', $U_ChiSquare);
        }
    }
    $workbook->close();
    print "$outxls\n";
}
################################################################################

1;

################################################################################

sub count_prefunpref {
    my ($seq1, $seq2, $puntable, $CODONS, $synchanges) = @_;
    my @puntable = split '', $puntable;
    my %results;
    my $gap_cnt = 0;

    my %input;
    my $seqlen = length($seq1);
    for (my $j=0; $j< $seqlen; $j+=3) {
#         ########################################
#         if ($j == 270000) {
#             print STDERR "STDERR pref/unpref from $seq1id to $seq2id\n";
#             for my $key (sort keys %results) {
#                 print STDERR "STDERR $key = $results{$key}\n";
#             }
#         }
#         ########################################
	$input{'cod1'} = substr($seq1, $j,3);
	$input{'cod2'} = substr($seq2, $j,3);
        $results{'00_total'}++;
	#ignore codon if being compared with gaps! 
	if ($input{'cod1'} =~ /\-/ || $input{'cod2'} =~ /\-/) {
	    $gap_cnt += 3; #just increments once if there is a pair of gaps
	    next;
	}

	my ($diff_cnt, $same) = count_diffs(\%input);

	#ignore if codons are identical
	if ($diff_cnt == 0) {
            $results{'01_no_nt_change'}++;
        } elsif ($diff_cnt == 1) {
        #get pun values for each codon
            if (1 == $synchanges->{$input{'cod1'}}{$input{'cod2'}}) {
                my $cod1 = $input{'cod1'};
                my $cod2 = $input{'cod2'};
                my $pun_cod1 = $puntable[$CODONS->{$cod1}];
                my $pun_cod2 = $puntable[$CODONS->{$cod2}];
                my $pun_string = "11_" . $pun_cod1 . $pun_cod2;
                if (uc($pun_cod1) eq uc($pun_cod2)) {
                    $results{'08_no_state_change'}++;
                } else {
                    $results{'09_state_change'}++;
                };
                $results{$pun_string}++;
                $results{'07_single_syn'}++;
                if (($pun_cod2 =~ /p/i) && (!($pun_cod1 =~ /p/i))) {
                    $results{'10-p'}++;
                } elsif (($pun_cod2 =~ /u/i) && (!($pun_cod1 =~ /u/i))) {
                    $results{'10-u'}++;
                } elsif (($pun_cod2 =~ /n/i) && (!($pun_cod1 =~ /n/i))) {
                    $results{'10-n'}++;
                } elsif (($pun_cod2 =~ /a/i) && (!($pun_cod1 =~ /a/i))) {
                    $results{'10-a'}++;
                } else {
                    my $string = "10__" . $pun_string;
                    $results{$pun_string}++;
                }
            } else {
                #ignore if codons are nonsyn
                $results{'06_single_nsyn'}++;
            }
        } elsif ($diff_cnt == 2) {
                $results{'04_double'}++;
        } elsif ($diff_cnt == 3) {
                $results{'03_triple'}++;
        } else {
            warn("Unrecognized type of change\n");
        }
    }
    $results{'05_single'} = $results{'07_single_syn'} + $results{'06_single_nsyn'};
    $results{'02_nt_change'} = ($results{'00_total'} - $results{'01_no_nt_change'});
    $results{'00_gaps'} = $gap_cnt;
    return %results;
}

################################################################################

sub count_diffs {
    #counts the number of nucleotide differences between 2 codons
    # returns this value plus the codon index of which nucleotide is the same when 2
    #nucleotides are different. This is so analyse_mutations() knows which nucleotides
    # to change.
    my $ref = shift;
    my $cnt = 0;
    my $same= undef;
    #just for 2 differences
    for (0..2) {
	if (substr($ref->{'cod1'}, $_,1) ne substr($ref->{'cod2'}, $_, 1)){
	    $cnt++;
	} else {
	    $same = $_;
	}
    }
    return ($cnt, $same);
}

################################################################################

1;
