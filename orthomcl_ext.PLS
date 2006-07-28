#!/usr/local/bin/perl

# POD documentation - main docs before the code

# Modified by Albert Vilella
# based on ORTHOMCL [2005-03-14] Version 1.2
# by Li Li, Feng Chen <fengchen@sas.upenn.edu>
#
# Copyright (C) 2004 by University of Pennsylvania, Philadelphia, PA
# USA.  All rights reserved.

=head1 NAME

orthomcl_ext.PLS - DESCRIPTION 

=head1 SYNOPSIS

There should be a BLOSUM62 file in the place where this script is called from.

perl orthomcl_ext.PLS --mode 1 -mcldir /home/avb/9_opl/mcl/mcl-04-314/src/shmcl/ -blastdir /home/avb/9_opl/blast/bin/ -datadir /home/avb/9_opl/orthomcl/sample_data -workdir /home/avb/ --fa_files Ath.fa,Sce.fa,Hsa.fa

perl orthomcl_ext.PLS --mode 2 --former_run_dir 20050908_230614 -mcldir /home/avb/9_opl/mcl/mcl-04-314/src/shmcl/ --inflation 1.5

perl orthomcl_ext.PLS --mode 3 --usr_blast_file /home/avb/20050908_230614/tmp/all.blast --usr_gg_file /home/avb/20050908_230614/tmp/all.gg -mcldir /home/avb/9_opl/mcl/mcl-04-314/src/shmcl/

perl orthomcl_ext.PLS --mode 4 --usr_bpo_file /home/avb/20050908_230614/tmp/all.bpo --usr_gg_file /home/avb/20050908_230614/tmp/all.gg -mcldir /home/avb/9_opl/mcl/mcl-04-314/src/shmcl/

perl orthomcl_ext.PLS --mode 5 --former_run_dir 20050908_230614 --usr_Taxa_file /home/avb/20050908_230614/tmp/all.gg --inflation=1.1 -mcldir /home/avb/9_opl/mcl/mcl-04-314/src/shmcl/

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
use File::Basename;
use Bio::SearchIO;
use Storable;
# use orthomcl_module; # All is contained here, no extra files

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $start_time = sprintf ("%04d%02d%02d_%02d%02d%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec);

my ($mode,$fa_files,$pv_cutoff,$pi_cutoff,
    $pmatch_cutoff,$inflation,$maximum_weight);
my ($usr_blast_file,$usr_bpo_file,$usr_gg_file,
    $usr_taxa_file,$former_run_dir,$email);    # For Mode 2, 3 or 4

my $command=basename($0)." ".join(' ',@ARGV)."\n";

my $BLASTALL                            = "/usr/local/bin/blastall";
my $FORMATDB                            = "/usr/local/bin/formatdb";
my $MCL                                 = "/usr/local/bin/mcl";
# path
my $PATH_TO_ORTHOMCL                    = "/home/avb/";   # must end with "/" 
my $ORTHOMCL_DATA_DIR                   = $PATH_TO_ORTHOMCL."/sample_data/";
my $ORTHOMCL_WORKING_DIR                = '';                            # will be changed with each run
my $ORTHOMCL_FORMER_RUN_DIR             = '';                            # if and only if user wants to use former run data
my $ORTHOMCL_TMP_DIR                    = '';                            # will be changed with each run

# cut-off variables
my $BLAST_PVALUE_CUTOFF_DEFAULT         = 1e-5;
my $PERCENT_IDENTITY_CUTOFF_DEFAULT     = 0;
my $PERCENT_MATCH_CUTOFF_DEFAULT        = 0;
my $MAX_WEIGHT_DEFAULT                  = 100;     # You need to see what's the second minimum p-value
                                                   # of your blast software. For example, if 9e-99 is 
                                                   # the case you should use -log(9e-99)=100.
my $MCL_INFLATION_DEFAULT               = 1.5;

# files
my $all_fa_file                         = "all.fa";
my $genome_gene_file                    = "all.gg";
my $blast_file                          = "all.blast";
my $bpo_file                            = "all.bpo";
my $bpo_idx_file                        = "all_bpo.idx";
my $bpo_se_file                         = "all_bpo.se";
my $matrix_file                         = "all_ortho.mtx";
my $index_file                          = "all_ortho.idx";
my $mcl_file                            = "all_ortho.mcl";
my $bbh_file                            = "all_blast.bbh";
my $mcl_bi_file                         = "all_orthomcl.out";
my $parameter_log_file                  = "parameter.log";
my $orthomcl_log_file                   = "orthomcl.log";

# global variables storing data
my @taxa=();
my %gindex=();     # taxon_id -> reference to the array with gene_id as elements
my %gindex2=();    # gene_id -> taxon_id
my %blastquery=();

my %graph=();
my %weight=();

my $address_ref='';
my $blastquery_ref='';
my @mcl_index=();

my ($blast_dir,$mcl_dir,$orthomcl_dir,$data_dir,$blast_cutoff);

&GetOptions(
			"mode=i"              => \$mode,
			"fa_files=s"          => \$fa_files,
			"pv_cutoff=s"         => \$pv_cutoff,
			"pi_cutoff=f"         => \$pi_cutoff,
			"pmatch_cutoff=f"     => \$pmatch_cutoff,
			"inflation=f"         => \$inflation,
			"maximum_weight=i"    => \$maximum_weight,
			"usr_blast_file=s"    => \$usr_blast_file,
			"usr_bpo_file=s"      => \$usr_bpo_file,
			"usr_gg_file=s"       => \$usr_gg_file,
			"usr_taxa_file=s"     => \$usr_taxa_file,
			"former_run_dir=s"    => \$former_run_dir,
			"blast|blastdir=s"    => \$blast_dir,
			"mcl|mcldir=s"        => \$mcl_dir,
			"work|workdir=s"      => \$orthomcl_dir,
			"data|datadir=s"      => \$data_dir,
			"blast_cutoff=s"      => \$blast_cutoff,
                        "email=s"             => \$email
);

$BLASTALL                            = "$blast_dir"."/blastall" if ($blast_dir);
$FORMATDB                            = "$blast_dir"."/formatdb" if ($blast_dir);
$MCL                                 = "$mcl_dir"."/mcl" if ($mcl_dir);
$ORTHOMCL_DATA_DIR                   = "$data_dir" if ($data_dir);
$PATH_TO_ORTHOMCL                    = "$orthomcl_dir" if ($orthomcl_dir);
$BLAST_PVALUE_CUTOFF_DEFAULT         = "$blast_cutoff" if ($blast_cutoff);

if (!defined $mode) {printHelp();}

#set the default
$pv_cutoff      = $pv_cutoff      ? $pv_cutoff      : $BLAST_PVALUE_CUTOFF_DEFAULT;
$pi_cutoff      = $pi_cutoff      ? $pi_cutoff      : $PERCENT_IDENTITY_CUTOFF_DEFAULT;
$pmatch_cutoff  = $pmatch_cutoff  ? $pmatch_cutoff  : $PERCENT_MATCH_CUTOFF_DEFAULT;
$inflation      = $inflation      ? $inflation      : $MCL_INFLATION_DEFAULT;
$maximum_weight = $maximum_weight ? $maximum_weight : $MAX_WEIGHT_DEFAULT;

my (%connect, %ortho);


if ($mode == 1) {
	if (defined $fa_files) {
		&constructDirectory($start_time);
		&constructAllFasta($fa_files,$all_fa_file);                                #construct all.fa file
		&executeFORMATDB($all_fa_file);                                            # and run blast
		&executeBLASTALL($all_fa_file,$blast_file,$all_fa_file,$pv_cutoff);
		&blast_parse($blast_file,$bpo_file,$pv_cutoff) unless (-e $bpo_file);
	} else {dieWithUnexpectedError("In Mode 1, NAMES OF FASTA FILES need to be given!");}
}
elsif ($mode == 2) {
	if (defined $former_run_dir) {
		&constructDirectory($start_time,$former_run_dir);
		&read_ggfile($genome_gene_file);
	} else {dieWithUnexpectedError("In Mode 2, FORMER RUN DIRECTORY needs to be given!");}
}
elsif ($mode == 3) {
	if ((defined $usr_blast_file) && (defined $usr_gg_file)) {
		&constructDirectory($start_time);
		$all_fa_file      = 'N/A';
		$genome_gene_file = $usr_gg_file;
		read_ggfile($genome_gene_file);
		$blast_file       = $usr_blast_file;
		&blast_parse($blast_file,$bpo_file,$pv_cutoff) unless (-e $bpo_file);
	} else {dieWithUnexpectedError("In Mode 3, BLAST OUT FILE and GENOME-GENE FILE are required!");}
}
elsif ($mode == 4) {
	if ((defined $usr_bpo_file) && (defined $usr_gg_file)) {
		&constructDirectory($start_time);
		$all_fa_file      = 'N/A';
		$genome_gene_file = $usr_gg_file;
		$blast_file       = 'N/A';
		read_ggfile($genome_gene_file);
		$bpo_file         = $usr_bpo_file;
		if ($usr_bpo_file =~ m/(\S+)\.(\S+)/) {
			$bpo_idx_file     = $1.'_bpo.idx';
			$bpo_se_file      = $1.'_bpo.se';
		} else {
			$bpo_idx_file     = $usr_bpo_file.'_bpo.idx';
			$bpo_se_file      = $usr_bpo_file.'_bpo.se';
		}
	} else {dieWithUnexpectedError("In Mode 4, BPO (BLAST PARSE OUT) FILE and GG (GENOME-GENE RELATION) FILE are required!");}
} 
elsif ($mode == 5) {
	if ((defined $former_run_dir) && (defined $usr_taxa_file)) {
		mode5($start_time,$command,$former_run_dir,$usr_taxa_file,$inflation);
                my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
                my $end_time = sprintf ("%04d%02d%02d_%02d%02d%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec);
		&write_endtime_in_parameter_log($end_time);
		write_log("\nStart Time: $start_time\nEnd Time:   $end_time\n");
		die "\nStart Time: $start_time\nEnd Time:   $end_time\n";

	} else {dieWithUnexpectedError("In Mode 5, FORMER RUN DIR and TAXA LIST FILE are required!");}
} 
else {dieWithUnexpectedError("Mode 1,2,3 or 4 needs to be given!");}


&write_parameter_log($start_time,$command,$mode,$pv_cutoff,$pi_cutoff,$pmatch_cutoff,$inflation,$maximum_weight);

&constructIDX_for_bpofile($bpo_file,$bpo_idx_file) unless (-e $bpo_idx_file);
&constructSE_for_bpofile($bpo_file,$bpo_se_file) unless (-e $bpo_se_file);
&open_bpofile($bpo_file);
&retrieve_from_file($bpo_idx_file,$bpo_se_file);


foreach my $taxon (@taxa) {
	print STDERR "\nIdentifying inparalogs from $taxon\n";
	write_log("\nIdentifying inparalogs from $taxon\n");
	@{$connect{$taxon.' '.$taxon}}  = &makeInparalog($taxon);                      # identification of inparalogs
}

for(my $i=0;$i<scalar(@taxa)-1;$i++) {
	for(my $j=$i+1;$j<scalar(@taxa);$j++) {
		print STDERR "\nIdentifying ortholog pairs between $taxa[$i] and $taxa[$j]\n";  
		write_log("\nIdentifying ortholog pairs between $taxa[$i] and $taxa[$j]\n");  
		@{$connect{$taxa[$i].' '.$taxa[$j]}} = &makeOrtholog($taxa[$i],$taxa[$j]); # identification of orthologs
		my %e = %{$connect{$taxa[$i].' '.$taxa[$j]}->[0]};
		my %w =  %{$connect{$taxa[$i].' '.$taxa[$j]}->[1]};
		my $sumw =  $connect{$taxa[$i].' '.$taxa[$j]}->[3];
		my $c = $connect{$taxa[$i].' '.$taxa[$j]}->[4];
		my %p1 = %{$connect{$taxa[$i].' '.$taxa[$i]}->[0]};
		my %p2 = %{$connect{$taxa[$j].' '.$taxa[$j]}->[0]};
		my %para;
		foreach my $p (keys %{$connect{$taxa[$i].' '.$taxa[$i]}->[0]}) {
			push (@{$para{$p}}, @{$p1{$p}}); }
		foreach my $p (keys %{$connect{$taxa[$j].' '.$taxa[$j]}->[0]}) {
			push (@{$para{$p}}, @{$p2{$p}}); }
		
		foreach my $n (keys %e) {
			$ortho{$n} = 1;
			my (@nodes1, @nodes2);

			if (exists($para{$n})) {push (@nodes1, $n, @{$para{$n}});}
			else {push (@nodes1, $n);}

			foreach (@{$e{$n}}) {
				if (exists($para{$_})) {push (@nodes2, $_, @{$para{$_}});}
				else {push (@nodes2, $_);}
			}

			@nodes1=@{nonredundant_list(\@nodes1)}; #can be commented
			@nodes2=@{nonredundant_list(\@nodes2)};

			for(my $k=0;$k<scalar(@nodes1);$k++) {
				for(my $l=0;$l<scalar(@nodes2);$l++) {
					
					next if(exists($w{$nodes1[$k].' '.$nodes2[$l]}));
					my ($pv1, $pv2);

					if (blastqueryab($nodes1[$k],$nodes2[$l])) {
						my ($s,$pm,$pe,$pi)=(blastqueryab($nodes1[$k],$nodes2[$l]))[0,3,4,5];
						next if($pm.'e'.$pe > $pv_cutoff || $pi< $pi_cutoff);
						if($pmatch_cutoff) {
							next if(&simspan($s) < $pmatch_cutoff);
						}
						if($pm==0) { $pv1 = $maximum_weight;} else { $pv1 = -log($pm.'e'.$pe)/log(10); }
					} else {next;}

					if (blastqueryab($nodes2[$l],$nodes1[$k])) {
						my ($s,$pm,$pe,$pi)=(blastqueryab($nodes2[$l],$nodes1[$k]))[0,3,4,5];
						next if($pm.'e'.$pe > $pv_cutoff || $pi < $pi_cutoff);
						if($pmatch_cutoff) {
							next if(&simspan($s) < $pmatch_cutoff);
						}
						if($pm==0) { $pv2 = $maximum_weight;} else { $pv2 = -log($pm.'e'.$pe)/log(10); }
						my $edge_ref=$connect{$taxa[$i].' '.$taxa[$j]}->[0];
						push (@{$edge_ref->{$nodes1[$k]}}, $nodes2[$l]);
						push (@{$edge_ref->{$nodes2[$l]}}, $nodes1[$k]);
						my $wt = ($pv1+$pv2)/2;
						# use averaged score as edge weight
						$w{$nodes1[$k].' '.$nodes2[$l]} = sprintf("%.3f", $wt);
						$w{$nodes2[$l].' '.$nodes1[$k]} = sprintf("%.3f", $wt);
						$sumw += $wt;
						$c++;
					}
				}
			}
		}
		my $avgw;
                if ($sumw <= 0) {
                    $avgw = 0;
                } else {
                    $avgw = $sumw/$c;
                }
		print STDERR "$taxa[$i] and $taxa[$j] average weight: $avgw\n";
		write_log("$taxa[$i] and $taxa[$j] average weight: $avgw\n");
		foreach my $p (keys %w) {
			$w{$p} = sprintf("%.3f", $w{$p}/$avgw);
		}
		$connect{$taxa[$i].' '.$taxa[$j]}->[1] = \%w;
	}
}

%blastquery=();
%gindex=();

foreach my $taxon (@taxa) {
	print STDERR "\ncalculate average weight from $taxon\n";
	write_log("\ncalculate average weight from $taxon\n");
	my %e = %{$connect{$taxon.' '.$taxon}->[0]};
	my %w = %{$connect{$taxon.' '.$taxon}->[1]};

	my $count=0; my $sum=0;
	foreach my $pair (keys %w) {
		my ($n,$p) = split(' ',$pair);
		next unless($ortho{$n} || $ortho{$p});
		$count++;
		$sum += $w{$n.' '.$p};
	}
        my $avgw;
	if ($sum <= 0) {
            $avgw = 0;
        } else {
             $avgw = $sum/$count;
        }

	print STDERR "$taxon average weight: $avgw\n";
	write_log("$taxon average weight: $avgw\n");
	foreach my $p (keys %w) {
            if ($avgw <= 0) {
                $w{$p} = sprintf("%.3f", "0");
            } else {
		$w{$p} = sprintf("%.3f", $w{$p}/$avgw);
            }
	}
	$connect{$taxon.' '.$taxon}->[1] = \%w; 
}
%ortho=();


foreach my $p (keys %connect) {
	my %e = %{$connect{$p}->[0]};
	my %w =  %{$connect{$p}->[1]};
	
	foreach my $n (keys %e) {
		push(@{$graph{$n}}, @{$e{$n}});
		delete $e{$n};
	}
	%e=();
	foreach my $n (keys %w) {
		$weight{$n} = $w{$n};
		delete $w{$n};
	}
	%w=();
	delete $connect{$p};
}
%connect=();

write_matrix_index($matrix_file,$index_file);

%graph=();
%weight=();

executeMCL($matrix_file,$mcl_file,$inflation);
mcl_backindex($mcl_file,$mcl_bi_file);
%gindex2=();

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $end_time = sprintf ("%04d%02d%02d_%02d%02d%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec);
print STDERR "\nStart Time: $start_time\nEnd Time:   $end_time\n";
write_log("\nStart Time: $start_time\nEnd Time:   $end_time\n");
&write_endtime_in_parameter_log($end_time);

if ($email) {
    eval 'use Mail::Sendmail';
    if ($@) {
        print "Mail::Sendmail required for email sending";
    } else {
        my $message = "Start Time: $start_time\nEnd Time:   $end_time\n";
        my %mail = ( To      => "$email",
                     From    => 'orthomcl@localhost',
                     Subject => "Sendmail batch: $0 job finished.",
                     Message => "$message"
                   );

        sendmail(%mail) or die $Mail::Sendmail::error;
        print STDERR "Mail sent. Log says:\n", $Mail::Sendmail::log;
        print STDERR "\n";
    }
}


#######################################SUBROUTINES###########################################
# This subroutine is an important part of OrthoMCL, used to
# look for inparalog (recent paralog) which is defined as 
# reciprocal better hits here.
# Please refer to the OrthoMCL paper for more details.
# One Arguments:
# 1. String Variable: Taxon name
# Last modified: 07/19/04
sub makeInparalog {
	my $taxon = $_[0];
	my (%seqs, %inbest, %pvalue,%sim);
	foreach (@{$gindex{$taxon}}) {$seqs{$_} = 1;}
	foreach my $qid (keys %seqs) {
		my ($sStart,$sEnd);
		if (defined $blastquery{$qid}) {
			($sStart,$sEnd)=split (";",$blastquery{$qid});
		} else {next;}
		my @sorted_simid=pvtie_sort($sStart,$sEnd,$taxon);
		LINE:foreach (0..$#sorted_simid) {
			my ($s,$sid,$pm,$pe,$pi)=(&getline_from_bpofile($sorted_simid[$_]))[0,3,5,6,7];
			if ($sid ne $qid) {
				last LINE unless ($seqs{$sid});                                  ## better hit not from the same species
				last LINE if($pm.'e'.$pe > $pv_cutoff || $pi < $pi_cutoff);      ## better hit not meet the cutoff
				if($pmatch_cutoff) {
					next LINE if(&simspan($s) < $pmatch_cutoff);
				}
				push(@{$inbest{$qid}}, $sid);
				$pvalue{$qid.' '.$sid} = $pm.'e'.$pe; 
			}
		}
	}
	my @b = keys %inbest;
	print STDERR scalar(@b)." sequences have better hits within species\n";
	write_log(scalar(@b)." sequences have better hits within species\n");
	return &matrix(\%inbest, \%pvalue);
} ##makeInparalog




# This subroutine is an important part of OrthoMCL, used to
# look for ortholog which is defined as the reciprocal best
# hit between two species.
# Please refer to the OrthoMCL paper for more details.
# Two Arguments:
# 1. String Variable: Taxon name
# 2. String Variable: Taxon name
# Last modified: 07/19/04
sub makeOrtholog {
	my ($ta,$tb) = @_;
	my (@seqs,%best,%sim,%pvalue);
	foreach my $qid (@{$gindex{$ta}}) {
		my ($sStart,$sEnd);
		if (defined $blastquery{$qid}) {
			($sStart,$sEnd)=split (";",$blastquery{$qid});
		} else {next;}
		my ($lastpm,$lastpe);
		my $hit_id=0;
		LINE:foreach ($sStart..$sEnd) {
			my ($s,$sid,$pm,$pe,$pi)=(&getline_from_bpofile($_))[0,3,5,6,7];
			if (defined $gindex2{$sid}) {
				if ($gindex2{$sid} eq $tb) {
					$hit_id++;
					if ($hit_id==1) {
						push(@{$sim{$qid}},"$sid,$pm,$pe,$pi,$s");
						$lastpm=$pm;$lastpe=$pe;
					}
					else {
						if (($lastpm==$pm) && ($lastpe==$pe)) {
							push(@{$sim{$qid}},"$sid,$pm,$pe,$pi,$s");
						}
						else {last LINE;}
					}
				}
			} else {print STDERR "$sid gindex2 not defined; lineid: $_\n";}
		}
	}
	foreach my $qid (@{$gindex{$tb}}) {
		my ($sStart,$sEnd);
		if (defined $blastquery{$qid}) {
			($sStart,$sEnd)=split (";",$blastquery{$qid});
		} else {next;}
		my ($lastpm,$lastpe);
		my $hit_id=0;
		LINE:foreach ($sStart..$sEnd) {
			my ($s,$sid,$pm,$pe,$pi)=(&getline_from_bpofile($_))[0,3,5,6,7];
			if (defined $gindex2{$sid}) {
				if ($gindex2{$sid} eq $ta) {
					$hit_id++;
					if ($hit_id==1) {
						push(@{$sim{$qid}},"$sid,$pm,$pe,$pi,$s");
						$lastpm=$pm;$lastpe=$pe;
					}
					else {
						if (($lastpm==$pm) && ($lastpe==$pe)) {
							push(@{$sim{$qid}},"$sid,$pm,$pe,$pi,$s");
						}
						else {last LINE;}
					}
				}
			} else {print STDERR "$sid gindex2 not defined; lineid: $_\n";}
		}
	}

	foreach my $q (keys %sim) {
		foreach (@{$sim{$q}}) {
			my @bla=split (',',$_);
			next if($bla[1].'e'.$bla[2] > $pv_cutoff || $bla[3]< $pi_cutoff);
			if($pmatch_cutoff) {
				next if(&simspan($bla[4]) < $pmatch_cutoff);
			}
			push(@{$best{$q}}, $bla[0]);
			$pvalue{$q.' '.$bla[0]} = $bla[1].'e'.$bla[2];

		}
	}
	my @b = keys %best;
	print STDERR scalar(@b)." sequences have best hits from the other species\n";
	write_log(scalar(@b)." sequences have best hits from the other species\n");
	return &matrix(\%best, \%pvalue);
} ## makeOrtholog




# This subroutine is used to choose two-way hits among one-way hits (best
# hits between two species or better hits within one species), 
# calculate the weight between two nodes (minus logrithm of the p-value, 
# or $MAX_WEIGHT_DEFAULT for p-value 0 ), and calculate average
# weight among all inparalogs within one species or all orthologs between
# two species. (Weighting process takes place in the main script)
# Two Arguments:
# 1. Reference Variable: reference to a hash which stores all the possible
#    gene pairs (one-way best hit, or better hit).
# 2. Reference Variable: reference to a hash which stores the pvalue for
#    the gene pairs.
# Last modified: 07/20/04
sub matrix {
	my %best      = %{$_[0]};
	my %pvalue    = %{$_[1]};
	my (%edge, %weight);
	my $count=0;
	my $sumw=0;

	foreach my $query (sort keys %best) {
		foreach my $subject (@{$best{$query}}) {
			next if($weight{$query.' '.$subject});
			my $flag = 0;
			foreach my $q (@{$best{$subject}}) {
				if($q eq $query) { $flag = 1; }
			}
			if($flag == 1) {
				push (@{$edge{$query}}, $subject);
				push (@{$edge{$subject}}, $query);
				#use -logP as weights and treat P=0 as -logP=$maximum_weight (DEFAULT=300)
				my ($pv1, $pv2);
				if($pvalue{$query.' '.$subject} == 0) {
					$pv1 = $maximum_weight;
				}else { 
					$pv1 = -log($pvalue{$query.' '.$subject})/log(10);
				}	    
				if($pvalue{$subject.' '.$query} == 0) {
					$pv2 = $maximum_weight;
				}else {
					$pv2 = -log($pvalue{$subject.' '.$query})/log(10);
				}
				write_bbh("$query	$subject	".$pvalue{$query.' '.$subject}."	".$pvalue{$subject.' '.$query}."\n");
				my $w = ($pv1+$pv2)/2;
				$sumw += $w;
				$count++;
				# use averaged score as edge weight
				$weight{$query.' '.$subject} = sprintf("%.3f", $w);
				$weight{$subject.' '.$query} = sprintf("%.3f", $w);
			}
		}
	}
        my $avgw;
	if ($sumw <= 0) {
            $avgw = 0;
        } else {
             $avgw = $sumw/$count;
        }

	return (\%edge, \%weight, $avgw, $sumw, $count);
} ## matrix


# This subroutine is used by the subroutine makeInparalog,
# to solve the pv_tie problem. It rearranges the pv-tied blast hits
# so that the hits from a specific taxon are moved higher than hits from
# other species.
# Three Arguments:
# 1. Number Variable: starting line id (or similarity id) of bpo file (blast
#    parse out file)
# 2. Number Variable: ending line id (or similarity id) of bpo file (blast
#    parse out file)
# 3. String Variable: taxon
# Last modified: 07/20/04
sub pvtie_sort {
	my ($sStart,$sEnd,$taxon)=@_;
	my (@sorted_simid,@tmp);
	my ($lastpm,$lastpe)=(getline_from_bpofile($sStart))[5,6];
	foreach ($sStart..$sEnd) {
		my ($s,$sid,$pm,$pe,$pi)=(&getline_from_bpofile($_))[0,3,5,6,7];
		if (($lastpm==$pm) && ($lastpe==$pe)) {
			if ($gindex2{$sid} eq $taxon) {
				push (@sorted_simid,$s);
			}
			else {
				push (@tmp,$s);
			}
		}
		else {
			if (scalar(@tmp)>0) {
				push (@sorted_simid,@tmp);
				@tmp=();
			} 
			if ($gindex2{$sid} eq $taxon) {
				push (@sorted_simid,$s);
			}
			else {
				push (@tmp,$s);
			}
		}
		$lastpm=$pm;$lastpe=$pe;
	}
	if (scalar(@tmp)>0) {push (@sorted_simid,@tmp);}
	return @sorted_simid;
} ## pvtie_sort





# This subroutine, together with matchlen, are used to calculate
# how much of the query sequences match each other.
# One Argument:
# 1. Number Variable: line id (or similarity id) of bpo file (blast
# parse out file)
# Last modified: 07/21/04
sub simspan {
	my $s = $_[0];
	my (%sub_start, %sub_length, %query_start, %query_length);
	my @hsp=split ('\.',(&getline_from_bpofile($s))[8]);
	foreach (@hsp) {
		if (/(\d+)\:(\d+)\-(\d+)\:(\d+)\-(\d+)/) {
			$sub_start{$1}=$4; 
			$sub_length{$1}=$5-$4+1;
			$query_start{$1}=$2;
			$query_length{$1}=$3-$2+1;
		}
	}
	my $match_lengths = &matchlen(\%sub_start,\%sub_length);
	my $match_lengthq = &matchlen(\%query_start,\%query_length);			
	my ($lengthq,$lengths)=(&getline_from_bpofile($s))[2,4];   # June 3
	if($lengths >= $lengthq) {
		return 100*$match_lengthq/$lengthq;
	}else{
		return 100*$match_lengths/$lengths;
	}
} ##simspan





# This subroutine, together with simspan, are used to calculate
# how much of the query sequences match each other.
# Two Arguments:
# 1. Reference Variable: reference to an hash which stores the starting
#    position of each HSP.
# 2. Reference Variable: reference to an hash which stores the length
#    of each HSP.
# Last modified: 07/19/04
sub matchlen {

	my %start        = %{$_[0]}; 
	my %length       = %{$_[1]};

	my @starts = sort{$start{$a}<=>$start{$b}} (keys %start);
	return $length{$starts[0]} if(scalar(@starts)==1);
	my $i=1; 
	my  $match_length = $length{$starts[0]}; 
	my $pos = $length{$starts[0]} + $start{$starts[0]} ;
	while($i<scalar(@starts)) {

	if($length{$starts[$i]} + $start{$starts[$i]} <= $pos) {
		$i++;
		next;
	}
	if($start{$starts[$i]}> $pos) {
		$match_length += $length{$starts[$i]};
		$pos = $start{$starts[$i]} + $length{$starts[$i]};
	}else {
		$match_length += $length{$starts[$i]} - ($pos - $start{$starts[$i]});
		$pos = $start{$starts[$i]} + $length{$starts[$i]};
	}
	$i++;
	}

	return $match_length;
} ## matchlen



# Last modified: 07/22/04
sub printHelp {
	my (@foo) = <DATA>;
	print @foo;
	exit 1;
}

# softwares

# This subroutine is used to construct the directories to store
# the intermediate files and the final files.
# One or Two arguments:
# 1. String Variable: `date`
# 2. String Variable: (Optional) the existing directory, e.g. "July_21_2" where blast
#    file and bpo files are located. This option is useful with species
#    list unchanged (blast are already performed) but with orthomcl parameters
#    being changed.
# Last modified: 07/21/04
sub constructDirectory {
	my $date                 = $_[0];
	my $former_run_dir       = $_[1];

	$ORTHOMCL_WORKING_DIR = $PATH_TO_ORTHOMCL."$date";
	my $no=1;$no = sprintf("%03d",$no);
	if (-e $ORTHOMCL_WORKING_DIR) {
		$no++;$no = sprintf("%03d",$no);
		while (-e $ORTHOMCL_WORKING_DIR."_".$no) {
			$no++;$no = sprintf("%03d",$no);
		}
		$ORTHOMCL_WORKING_DIR.="_$no";
	}
	$ORTHOMCL_WORKING_DIR.="/";
	system ("mkdir $ORTHOMCL_WORKING_DIR");
	$mcl_bi_file=$ORTHOMCL_WORKING_DIR.$mcl_bi_file;
	$parameter_log_file=$ORTHOMCL_WORKING_DIR.$parameter_log_file;
# 	$orthomcl_log_file=$ORTHOMCL_WORKING_DIR.$orthomcl_log_file;open (LOG,">$orthomcl_log_file") or die "can't open file";
	open (LOG,">$orthomcl_log_file") or die "can't open file";
# 	$bbh_file=$ORTHOMCL_WORKING_DIR.$bbh_file;open (BBH,">$bbh_file") or die "can't create $bbh_file";
	open (BBH,">$bbh_file") or die "can't create $bbh_file";
	print STDERR "### WORKING DIRECTORY: ###\n  $ORTHOMCL_WORKING_DIR\n\n";
#	print LOG "### WORKING DIRECTORY: ###\n  $ORTHOMCL_WORKING_DIR\n\n";
	write_log("### WORKING DIRECTORY: ###\n  $ORTHOMCL_WORKING_DIR\n\n");

	$ORTHOMCL_TMP_DIR=$ORTHOMCL_WORKING_DIR."tmp/";
	system ("mkdir $ORTHOMCL_TMP_DIR");
	
	if (defined $former_run_dir) {
		$ORTHOMCL_FORMER_RUN_DIR = $PATH_TO_ORTHOMCL.$former_run_dir."/tmp/";
		$all_fa_file             = $ORTHOMCL_FORMER_RUN_DIR.$all_fa_file;
		$genome_gene_file        = $ORTHOMCL_FORMER_RUN_DIR.$genome_gene_file;
		$blast_file              = $ORTHOMCL_FORMER_RUN_DIR.$blast_file;
		$bpo_file                = $ORTHOMCL_FORMER_RUN_DIR.$bpo_file;
		$bpo_idx_file            = $ORTHOMCL_FORMER_RUN_DIR.$bpo_idx_file;
		$bpo_se_file             = $ORTHOMCL_FORMER_RUN_DIR.$bpo_se_file;
	} else {
		$all_fa_file             = $ORTHOMCL_TMP_DIR.$all_fa_file;
		$genome_gene_file        = $ORTHOMCL_TMP_DIR.$genome_gene_file;
		$blast_file              = $ORTHOMCL_TMP_DIR.$blast_file;
		$bpo_file                = $ORTHOMCL_TMP_DIR.$bpo_file;
		$bpo_idx_file            = $ORTHOMCL_TMP_DIR.$bpo_idx_file;
		$bpo_se_file             = $ORTHOMCL_TMP_DIR.$bpo_se_file;
	}
	$matrix_file       = $ORTHOMCL_TMP_DIR.$matrix_file;
	$index_file        = $ORTHOMCL_TMP_DIR.$index_file;
	$mcl_file          = $ORTHOMCL_TMP_DIR.$mcl_file;
}

sub write_log {
	my $comment = $_[0];
	print LOG $comment;
}
sub write_bbh {
	my $comment = $_[0];
	print BBH $comment;
}


# This subroutine is used to construct the directories to store
# the intermediate files and the final files.
# One or Two arguments:
# 1. String Variable: `date`
# 2. String Variable: the existing directory, e.g. "July_21_2" where matrix
#    file and index file are located. This option is useful with species
#    list changed but other orthomcl parameters unchanged.
# Last modified: 11/22/04
sub mode5 {
	my $date                 = $_[0];
	my $command              = $_[1];
	my $former_run_dir       = $_[2];
	my $usr_taxa_file        = $_[3];
	my $inflation            = $_[4];


	$ORTHOMCL_WORKING_DIR = $PATH_TO_ORTHOMCL."$date";
	my $no=1;$no = sprintf("%03d",$no);
	if (-e $ORTHOMCL_WORKING_DIR) {
		$no++;$no = sprintf("%03d",$no);
		while (-e $ORTHOMCL_WORKING_DIR."_".$no) {
			$no++;$no = sprintf("%03d",$no);
		}
		$ORTHOMCL_WORKING_DIR.="_$no";
	}
	$ORTHOMCL_WORKING_DIR.="/";
	system ("mkdir $ORTHOMCL_WORKING_DIR");

	$ORTHOMCL_TMP_DIR=$ORTHOMCL_WORKING_DIR."tmp/";
	system ("mkdir $ORTHOMCL_TMP_DIR");
	
	my ($former_matrix_file,$former_index_file);
	if (defined $former_run_dir) {
		$ORTHOMCL_FORMER_RUN_DIR = $PATH_TO_ORTHOMCL.$former_run_dir.'/';

		open(FORMERLOG,$ORTHOMCL_FORMER_RUN_DIR.$parameter_log_file) || die;
		while (<FORMERLOG>) {
			if (/all\_fa\_file\s+(\S+)/) {
				my $file=$1;
				if ($file=~/\//) {$all_fa_file=$file;}
				else {$all_fa_file=$PATH_TO_ORTHOMCL.$file;}
			} elsif (/genome\_gene\_file\s+(\S+)/) {
				my $file=$1;
				if ($file=~/\//) {$genome_gene_file=$file;}
				else {$genome_gene_file=$PATH_TO_ORTHOMCL.$file;}
			} elsif (/blast\_file\s+(\S+)/) {
				my $file=$1;
				if ($file=~/\//) {$blast_file=$file;}
				else {$blast_file=$PATH_TO_ORTHOMCL.$file;}
			} elsif (/blast\_parse\_out\_file\s+(\S+)/) {
				my $file=$1;
				if ($file=~/\//) {$bpo_file=$file;}
				else {$bpo_file=$PATH_TO_ORTHOMCL.$file;}
			} elsif (/bpo\_idx\_file\s+(\S+)/) {
				my $file=$1;
				if ($file=~/\//) {$bpo_idx_file=$file;}
				else {$bpo_idx_file=$PATH_TO_ORTHOMCL.$file;}
			} elsif (/bpo\_se\_file\s+(\S+)/) {
				my $file=$1;
				if ($file=~/\//) {$bpo_se_file=$file;}
				else {$bpo_se_file=$PATH_TO_ORTHOMCL.$file;}
			}
		}
		$bbh_file                = $ORTHOMCL_FORMER_RUN_DIR.$bbh_file;
		$ORTHOMCL_FORMER_RUN_DIR.='tmp/';
		$former_matrix_file      = $ORTHOMCL_FORMER_RUN_DIR.$matrix_file;
		$former_index_file       = $ORTHOMCL_FORMER_RUN_DIR.$index_file;
	} else {die "wrong subroutine call!\n";}
	$matrix_file       = $ORTHOMCL_TMP_DIR.$matrix_file;
	$index_file        = $ORTHOMCL_TMP_DIR.$index_file;
	$mcl_file          = $ORTHOMCL_TMP_DIR.$mcl_file;

	$mcl_bi_file=$ORTHOMCL_WORKING_DIR.$mcl_bi_file;
	$parameter_log_file=$ORTHOMCL_WORKING_DIR.$parameter_log_file;
	$orthomcl_log_file=$ORTHOMCL_WORKING_DIR.$orthomcl_log_file;open (LOG,">$orthomcl_log_file") or die "can't open file";

	print STDERR "### WORKING DIRECTORY: ###\n  $ORTHOMCL_WORKING_DIR\n\n";
	write_log("### WORKING DIRECTORY: ###\n  $ORTHOMCL_WORKING_DIR\n\n");

	read_ggfile($genome_gene_file);
	write_parameter_log($date,$command,5,"N/A","N/A","N/A",$inflation,"N/A");
	my %mcl_backindex;
	# read usr_taxa_file
	open (TAXA,$usr_taxa_file) or die "can't open $usr_taxa_file";
	while (<TAXA>) {
		if (/^([A-z]+)[\(\:\n]/) {
			my $taxon = $1;
			die "The taxon $1 is not defined!" unless (defined $gindex{$taxon});
			push(@mcl_index,@{$gindex{$taxon}});
		}
	}
	close(TAXA);
	for (my $i=0;$i<=$#mcl_index;$i++) {$mcl_backindex{$mcl_index[$i]} = $i;}

	# read former run index
	my @former_mcl_index;
	open (IDX,$former_index_file) or die "can't open $former_index_file";
	while (<IDX>) {
		if (/(\d+)	(\S+)/) {$former_mcl_index[$1]=$2;}
	}
	close(IDX);
	# read former run matrix
	my @matrix;
	open (MTX,$former_matrix_file) or die "can't open $former_matrix_file";
	while (<MTX>) {
		if (/^(\d+)\s+(.*)\$$/) {
			my $idx = $1;
			next unless (defined $mcl_backindex{$former_mcl_index[$idx]});
			my @score = split(" ",$2);
			foreach (@score) {
				if (/(\d+)\:(\S+)/) {
					next unless (defined $mcl_backindex{$former_mcl_index[$1]});
					push(@{$matrix[$mcl_backindex{$former_mcl_index[$idx]}]},$mcl_backindex{$former_mcl_index[$1]}.":".$2);
				}
			}
		}
	}
	close(MTX);

	my $size = scalar(@matrix);
	print STDERR "\nThere are $size sequences to cluster\n";
	write_log("\nThere are $size sequences to cluster\n");
	open (MTX,">$matrix_file") or dieWithUnexpectedError("cannot write to file $matrix_file");
	print MTX "(mclheader\nmcltype matrix\ndimensions ".$size."x".$size."\n)\n\n(mclmatrix\nbegin\n\n";

	for (my $i = 0; $i <= $#matrix;$i++) {
		print MTX $i."    ";
		foreach (@{$matrix[$i]}) {
			print MTX $_." ";
		}
		print MTX "\$\n";
	}
	print MTX ")\n\n";
	close (MTX);

	print STDERR "Matrix($size X $size) file $matrix_file generated\n";
	write_log("Matrix($size X $size) file $matrix_file generated\n");

	open(IDX,">$index_file") or dieWithUnexpectedError("cannot write to file $index_file");
	for (my $i = 0;$i <= $#mcl_index; $i++) {
		print IDX "$i\t$mcl_index[$i]\n";
	}
	close(IDX);
	print STDERR "\nIndex file $index_file generated\n";
	write_log("\nIndex file $index_file generated\n");
	
	executeMCL($matrix_file,$mcl_file,$inflation);
	mcl_backindex($mcl_file,$mcl_bi_file);
}


# Two arguments:
# 1. String Variable: Fasta file names, each representing one species, separated by comma, e.g. "Bme.fa,Ccr.fa,Eco.fa"
# 2. String Variable: Name of All_Fasta file
# Last modified: 07/19/04
sub constructAllFasta {
	my $fa_files=$_[0];
	my $all_fa_file=$_[1];

	my @fafiles=split (",",$fa_files);
	my $totalgeneno=0;
	open (ALLFA,">$all_fa_file") or dieWithUnexpectedError("can't write to $all_fa_file!");
	foreach my $fafile (sort @fafiles) {
		open (FA,$ORTHOMCL_DATA_DIR."/".$fafile) or dieWithUnexpectedError("can't open $fafile!");
		my $taxon=(split (".fa",$fafile))[0];
		push (@taxa,$taxon);
		while (<FA>) {
			$_=~s/\r|\n//g;
			print ALLFA "$_\n";
			if (/\>(\S+)/) {push (@{$gindex{$taxon}},$1);}
		}
		close (FA);
		foreach (@{$gindex{$taxon}}) {$gindex2{$_}=$taxon;}
		$totalgeneno+=@{$gindex{$taxon}};
	}
	close (ALLFA);

	open (GG,">$genome_gene_file");
	foreach my $taxon (@taxa) {
		print GG "$taxon:";
		foreach (@{$gindex{$taxon}}) {
			print GG " $_";
		}
		print GG "\n";
	}
	close (GG);

	print STDERR "\nFASTA file <$all_fa_file> (".@fafiles." genomes, $totalgeneno sequences) generated!\n\n";
	write_log("\nFASTA file <$all_fa_file> (".@fafiles." genomes, $totalgeneno sequences) generated!\n\n");
} ## constructAllFasta




# Two arguments:
# 1. String Variable: genome-gene file
# Last modified: 07/20/04
sub read_ggfile {
	my $gg_file=$_[0];
	my $totalgeneno=0;
	open (GG,$gg_file) or dieWithUnexpectedError("can't open $gg_file!");
	while (<GG>) {
		$_=~s/\r|\n//g;
		if (/(\S+)\(\d+\):(.*)/) {
			my $taxon=$1;
			my @genes=split (" ",$2);
			push (@taxa,$taxon);
			push (@{$gindex{$taxon}},@genes);
			foreach (@genes) {$gindex2{$_}=$taxon;}
			$totalgeneno+=scalar(@{$gindex{$taxon}});
		}
		elsif (/(\S+):(.*)/) {
			my $taxon=$1;
			my @genes=split (" ",$2);
			push (@taxa,$taxon);
			push (@{$gindex{$taxon}},@genes);
			foreach (@genes) {$gindex2{$_}=$taxon;}
			$totalgeneno+=scalar(@{$gindex{$taxon}});
		}
	}
	print STDERR "\nThere are ".@taxa." genomes, $totalgeneno sequences in $gg_file !\n\n";
	write_log("\nThere are ".@taxa." genomes, $totalgeneno sequences in $gg_file !\n\n");
	close (GG);
} ## read_ggfile




# One Argument:
# 1. String Variable: fasta file name
# Last modified: 07/19/04
sub executeFORMATDB {
	my $in = $_[0];
	print STDERR "\nRunning FormatDB\n";
	write_log("\nRunning FormatDB\n");
	print STDERR "  Native FormatDB command:\n  $FORMATDB -i $in -p T -a F -o F\n\n";
	write_log("  Native FormatDB command:\n  $FORMATDB -i $in -p T -a F -o F\n\n");
	system ("$FORMATDB -i $in -p T -a F -o F");
} ## executeFORMATDB




# Four or Five Arguments:
# 1. String Variable: fasta file name
# 2. String Variable: blast out file name
# 3. String Variable: database name
# 4. String Variable: pvalue cutoff
# 5. (Optional) String Variable: "+number", parallel computing
# Last modified: 07/19/04
sub executeBLASTALL {
	my $in        = $_[0];
	my $out       = $_[1];
	my $db        = $_[2];
	my $pv_cutoff = $_[3];
	print STDERR "\nRunning all-against-all BLAST\n";
	write_log("\nRunning all-against-all BLAST\n");
	print STDERR "  Native BLASTALL command:\n  $BLASTALL -p blastp -i $in -d $db -e $pv_cutoff -o $out\n";
	write_log("  Native BLASTALL command:\n  $BLASTALL -p blastp -i $in -d $db -e $pv_cutoff -o $out\n");
	system ("$BLASTALL -p blastp -i $in -d $db -e $pv_cutoff -o $out");
} ## executeBLASTALL




# Three Arguments:
# 1. String Variable: matrix file name
# 2. String Variable: mcl file name
# 3. Number Variable: inflation parameter for MCL algorithm
# Last modified: 07/19/04
sub executeMCL {
	my $matrix_file = $_[0];
	my $mcl_file    = $_[1];
	my $inflation  = $_[2];
	print STDERR "\nRun MCL program\n  $MCL $matrix_file -I $inflation -o $mcl_file\n";
	write_log("\nRun MCL program\n  $MCL $matrix_file -I $inflation -o $mcl_file\n");
	system("$MCL $matrix_file -I $inflation -o $mcl_file");
	if (-e $mcl_file) {
		print STDERR "\nMCL result $mcl_file generated!\n";
		write_log("\nMCL result $mcl_file generated!\n");}
	else {dieWithUnexpectedError("$mcl_file failed to be generated!");}
} ## executeMCL



####################
####################
#################### $pv_cutoff
# Two Arguments:
# 1. String Variable: blast out file
# 2. String Variable: blast parse out file
# Last modified: 07/19/04
sub blast_parse {
	my $blastfile        = $_[0];
	my $parseoutfile     = $_[1];
	my $pv_cutoff        = $_[2];

	open (PARSEOUT,">$parseoutfile");
	print STDERR "\nParsing blast result!\n";
	write_log("\nParsing blast result!\n");
	my $searchio = Bio::SearchIO->new(-file   => $blastfile,
									  -format => 'blast') or dieWithUnexpectedError("Blast parsing failed!");
	my $similarityid=1;
	while (my $result = $searchio->next_result()) {
		my $queryid=$result->query_name;
		my $querylen=$result->query_length;
		while( my $hit = $result->next_hit ) {
			next unless numeric_pvalue($hit->significance) < $pv_cutoff;
			my $subjectid=$hit->name;
			my $subjectlen=$hit->length;
			my $pvalue=numeric_pvalue($hit->significance);
			my ($pm,$pe);
			if ($pvalue==0) {$pm=0;$pe=0;}
			elsif ($pvalue=~/(\d+)e(\-\d+)/) {$pm=$1;$pe=$2;}
			my $simspanid=1;
			my $simspan='';
			my (@percentidentity,@hsplength);
			while( my $hsp = $hit->next_hsp ) {
				my $querystart=$hsp->start('query');
				my $queryend=$hsp->end('query');
				my $subjectstart=$hsp->start('sbjct');
				my $subjectend=$hsp->end('sbjct');
				$percentidentity[$simspanid]=$hsp->percent_identity;
				$hsplength[$simspanid]=$hsp->length('hit');
				$simspan.="$simspanid:$querystart-$queryend:$subjectstart-$subjectend.";
				$simspanid++;
			}
			my $sum_identical=0;
			my $sum_length=0;
			for (my $i=1;$i<$simspanid;$i++) {
				$sum_identical+=$percentidentity[$i]*$hsplength[$i];
				$sum_length+=$hsplength[$i];
			}
			my $percentIdent;
                        if ($sum_identical <= 0) {
                            $percentIdent = 0;
                        } else {
                            $percentIdent =int($sum_identical/$sum_length);
                        }
			print PARSEOUT "$similarityid;$queryid;$querylen;$subjectid;$subjectlen;$pvalue;$percentIdent;$simspan\n";
			$similarityid++;
		}
	}
	print STDERR "Parsing blast file finished\n";
	write_log("Parsing blast file finished\n");
	close(PARSEOUT);
} ## blast_parse




# This module is used to open BPO (Blast Parse Out) file, 
# to provide an filehandle for later use in blast query process.
# One Argument:
# 1. String Variable: BPO file name
# Last modified: 07/21/04
sub open_bpofile {
	my $bpo_file = $_[0];
	open (BPOFILE,$bpo_file) || dieWithUnexpectedError("can't open $bpo_file file"); # getline_from_bpofile subroutine 
																					 # will use this file handle
} # open_bpofile




# Make pvalue numeric, used by subroutine blast_parse
# One Arguments:
# 1. String Variable: pvalue
# Last modified: 07/19/04
sub numeric_pvalue {
	my $p=$_[0];
	if ($p=~/^e-(\d+)/) {return "1e-".$1;}
	else {return $p}
} # numeric_pvalue

# This subroutine is used to facilitate reading blast information
# from blastparseout file, instead of reading all blast result into
# memory.
# It works by storing every line (of blastparseout file)'s address
# into an array, and then uses Storable module to save this array
# into a file for future reference.
# Two Arguments:
# 1. String Variable: blast parse out file name
# 2. String Variable: index file name, storing index information of
#    blast parse out file.
# Last modified: 07/19/04
sub constructIDX_for_bpofile {
	my ($file,$idxfile)=@_;
	my @address;
	open (FILE,$file) or dieWithUnexpectedError("can't open $file file");
	push (@address, tell (FILE));
	while (<FILE>) {push (@address, tell (FILE));}
	Storable::store \@address, $idxfile;
	close (FILE);
} ## constructIDX_for_bpofile




# This subroutine is also used to facilitate reading blast information
# from blastparseout file, instead of reading all blast result into
# memory.
# It works by storing Starting and Ending line_ids of specific query_gene_id
# into a hash, and then uses Storable module to save this hash into a file
# for future reference.
# Two Arguments:
# 1. String Variable: blast parse out file name
# 2. String Variable: file name storing the hash
# Last modified: 07/19/04
sub constructSE_for_bpofile {
	my ($bpo_file,$bpo_se_file)=@_;
	my $lastquery='';my $lastsimid=0;
	open (BPOUT,$bpo_file) or dieWithUnexpectedError("can't open $bpo_file");
	while (<BPOUT>) {
		$_=~s/\r|\n//g;
		next unless ($_ ne '');
		my @bpout=split (";",$_);
		if ($lastquery ne $bpout[1]) {
			$blastquery{$bpout[1]}=$bpout[0];
			$blastquery{$lastquery}.=";".$lastsimid if $lastquery;
		}
		$lastquery=$bpout[1];
		$lastsimid=$bpout[0];
	}
	$blastquery{$lastquery}.=";".$lastsimid;
	Storable::store \%blastquery,$bpo_se_file;
	close(BPOUT);
} ## constructSE_for_bpofile




# This subroutine is used to retrieve @address and %blastquery
# from blast_idx_file and blast_SE_file
# Two arguments:
# 1. String Variable: blast_idx_file
# 2. String Variable: blast_SE_file
# Last modified: 07/19/04
sub retrieve_from_file {
	my $blastidxfile  = $_[0];
	my $blastsefile   = $_[1];
	$address_ref=retrieve ($blastidxfile);
	$blastquery_ref=retrieve ($blastsefile);
	%blastquery=%{$blastquery_ref};
} ## retrieve_from_file




# This subroutine is used to retrieve blast information from bpo file
# (blast parse out file), given the line_id.
# One Argument:
# 1. Reference Variable: reference to the address array
# 2. Number Variable: line_id of bpo file
# Last modified: 09/08/04
sub getline_from_bpofile {
	my $lineid         = $_[0];
	seek (BPOFILE,$address_ref->[$lineid-1],0);
	my $line=<BPOFILE>; 
	$line=~s/\r|\n//g;
	chop $line;
	my @bpout=split (";",$line);
	my ($pm,$pe);
	if ($bpout[5]==0) {$pm=0;$pe=0;}
#	elsif ($bpout[5]=~/(\d+)e(\-\d+)/) {$pm=$1;$pe=$2;}
	elsif ($bpout[5]=~/(\S+)e(\-\S+)/) {$pm=$1;$pe=$2;}  #For WU-BLAST their p-value has the pattern \d\.\de\-\d+  OR 0.
	else {$pm=$bpout[5];$pe=0;}
	return ($bpout[0],$bpout[1],$bpout[2],$bpout[3],,$bpout[4],$pm,$pe,$bpout[6],$bpout[7]);
} ## getline_from_bpofile




# This subroutine is used to retrieve blast information from bpo file
# (blast parse out file), given the query gene id and the subject gene
# id. If such information doesn't exist, zero is returned.
# Four Arguments:
# 1. String Variable: query gene id
# 2. String Variable: subject gene id
# Last modified: 07/19/04
sub blastqueryab {
	my $a              = $_[0];
	my $b              = $_[1];
	my ($s,$e);
	if (defined $blastquery_ref->{$a}) {
		($s,$e)=split (";",$blastquery_ref->{$a});
	} else {return 0;}
    foreach my $i ($s..$e) {
		my @c=(&getline_from_bpofile($i))[0,1,3,5,6,7];
		if ($c[2] eq $b) {return(@c);}
	}
	return 0;
} ## blastqueryab





# This subroutine is used to make a nonredundant list.
# One Argument:
# 1. Reference Variable: reference to an array
# Last modified: 07/19/04
sub nonredundant_list {
	my $list_ref=$_[0];
	my %nr=();
	foreach (@{$list_ref}) {$nr{$_}=1;}
	my @nr=sort (keys %nr);
	return \@nr;
} ## nonredundant_list





# This subroutine is used to generate input file (matrix file)
# for MCL and index file.
# Four arguments:
# 1. String Variable: matrix file name
# 2. String Variable: index file name
# Last modified: 02/24/05
sub write_matrix_index {

	my $matrix_file  = $_[0];
	my $index_file   = $_[1];

	my $size = scalar(keys %graph);
	print STDERR "\nThere are $size sequences to cluster\n";
	write_log("\nThere are $size sequences to cluster\n");
	open (MTX,">$matrix_file") or dieWithUnexpectedError("cannot write to file $matrix_file");
	print MTX "(mclheader\nmcltype matrix\ndimensions ".$size."x".$size."\n)\n\n(mclmatrix\nbegin\n\n";

	my $i=0;
	my %mcl_index2;
	foreach my $p (sort keys %graph) {
		$mcl_index2{$p}=$i;$mcl_index[$i]=$p;
		$i++;
	}
	foreach my $p (sort keys %graph) {

		print MTX $mcl_index2{$p}."    ";
		foreach my $m (@{$graph{$p}}) {
			print MTX $mcl_index2{$m}.":".$weight{$p.' '.$m}." ";
		}
		print MTX "\$\n";
	}
	print MTX ")\n\n";

	close (MTX);

	print STDERR "Matrix($size X $size) file $matrix_file generated\n";
	write_log("Matrix($size X $size) file $matrix_file generated\n");

	#################################WRITE INDEX FILE######################################
	open(IDX,">$index_file") or dieWithUnexpectedError("cannot write to file $index_file");
	foreach my $id (sort { $mcl_index2{$a} <=> $mcl_index2{$b} } keys %mcl_index2) {
		print IDX "$mcl_index2{$id}\t$id\n";
	}
	close(IDX);
	print STDERR "\nIndex file $index_file generated\n";
	write_log("\nIndex file $index_file generated\n");
}


# This subroutine is used to back index all gene_ids present 
# in MCL output and generate the final result.
# Two arguments:
# 1. String Variable: mcl file name
# 2. String Variable: mcl back_index file name
# Last modified: 07/20/04
sub mcl_backindex {
	
	my $mcl_file      = $_[0];
	my $mcl_bi_file   = $_[1];

	open (MCL,$mcl_file) or dieWithUnexpectedError("can't open $mcl_file");
	my $last=0;
	my $lastline='';
	my @mcl;
	while (<MCL>) {
# 		chomp;chop;
		chomp;
		if (/^(\d+)\s+(.*)\$/) {
			$mcl[$last]=$lastline;
			$last=$1;$lastline=$2;
		}
# 		elsif (/^(\d+)\s+(.*)/) {
# 			$mcl[$last]=$lastline;
# 			$last=$1;$lastline=$2;
# 		}
		elsif (/^\s+/) {$lastline.=$_;}
	}
	$mcl[$last]=$lastline;
	close (MCL);


	open (MCLBI,">$mcl_bi_file") or dieWithUnexpectedError("can't write to $mcl_bi_file");
	my $orthomcl_cluster_id=1;
	foreach my $mcl_cluster_id (0..$last) {
		$mcl[$mcl_cluster_id]=~s/\s+/ /g;
		my @g=split (" ",$mcl[$mcl_cluster_id]);
		my @taxa=();
		foreach (@g) {
			my $taxon=$gindex2{$mcl_index[$_]};
			my $presence=0;
			foreach (@taxa) {
				if ($_ eq $taxon) {$presence=1;}
			}
			push (@taxa,$taxon) unless ($presence);
		}

		print MCLBI "ORTHOMCL".$mcl_cluster_id."(".scalar(@g)." genes,".scalar(@taxa)." taxa):	";
		foreach (@g) {
			print MCLBI " $mcl_index[$_]($gindex2{$mcl_index[$_]})";
		}
		print MCLBI "\n";
#No Species Cutoff
	}

	print STDERR "\n\nFinal ORTHOMCL Result: $mcl_bi_file generated!!!\n\n";
	write_log("\n\nFinal ORTHOMCL Result: $mcl_bi_file generated!!!\n\n");

} ## mcl_backindex




# This subroutine is used to print out the parameters for future reference
# Eight Arguments:
# 1. String Variable: Start Time
# 2. String Variable: Command Line
# 3. Number Variable: OrthoMCL Mode number
# 4. Number Variable: P-value Cutoff
# 5. Number Variable: Percent Identity Cutoff
# 6. Number Variable: Percent Match Cutoff
# 7. Number Variable: MCL Inflation Parameter
# 8. Number variable: Maximum Weight
# Last modified: 07/21/04
sub write_parameter_log {
	my $starttime        = $_[0];
	my $command          = $_[1];
	my $mode             = $_[2];
	my $pv_cutoff        = $_[3];
	my $pi_cutoff        = $_[4];
	my $pmatch_cutoff    = $_[5];
	my $inflation        = $_[6];
	my $maximum_weight   = $_[7];
	open (PLOG,">$parameter_log_file") or die "can't open file";
	print PLOG "########################COMMAND######################\n";
	print PLOG "$command\n";
	print PLOG "######################PARAMETERS#####################\n";
	printf PLOG "%25s  %-10s\n","OrthoMCL Mode",$mode;
	printf PLOG "%25s  %-10s\n","P-value Cut-off",$pv_cutoff;
	printf PLOG "%25s  %-10s\n","Percent Identity Cut-off",$pi_cutoff;
	printf PLOG "%25s  %-10s\n","Percent Match Cut-off",$pmatch_cutoff;
	printf PLOG "%25s  %-10s\n","MCL Inflation",$inflation;
	printf PLOG "%25s  %-10s\n","Maximum Weight",$maximum_weight;
	printf PLOG "#######################SPECIES########################\n";
	foreach my $taxon (@taxa) {
		printf PLOG "%25s  %-20s\n",$taxon,scalar(@{$gindex{$taxon}})." genes";
	}
	printf PLOG "########################FILES########################\n";
	printf PLOG "%25s  %-50s\n","all_fa_file",$all_fa_file;
	printf PLOG "%25s  %-50s\n","genome_gene_file",$genome_gene_file;
	printf PLOG "%25s  %-50s\n","blast_file",$blast_file;
	printf PLOG "%25s  %-50s\n","blast_parse_out_file",$bpo_file;
	printf PLOG "%25s  %-50s\n","bpo_idx_file",$bpo_idx_file;
	printf PLOG "%25s  %-50s\n","bpo_se_file",$bpo_se_file;
	printf PLOG "%25s  %-50s\n","matrix_file",$matrix_file;
	printf PLOG "%25s  %-50s\n","index_file",$index_file;
	printf PLOG "%25s  %-50s\n","mcl_file",$mcl_file;
	printf PLOG "%25s  %-50s\n","mcl_bi_file",$mcl_bi_file;
	printf PLOG "%25s  %-50s\n\n","parameter_log_file",$parameter_log_file;
	print PLOG "######################START TIME#####################\n";
	print PLOG "$starttime\n";
	close (PLOG);
} ## write_parameter_log


# This subroutine is used to print out the end time in parameter log file
# One Argument:
# 1. String Variable: End Time
# Last modified: 09/03/04
sub write_endtime_in_parameter_log {
	my $endtime        = $_[0];
	open (PLOG,">>$parameter_log_file") or die "can't open file";
	print PLOG "########################END TIME#####################\n";
	print PLOG "$endtime\n";
	close (PLOG);
} ## write_endtime_in_parameter_log

# Last modified: 07/19/04
sub dieWithUnexpectedError {
	my $text = $_[ 0 ];
	die( "\n\n$0:\nUnexpected error (should not have happened):\n$text\n$!\n\n" );
} ## dieWithUnexpectedError


1;


######################################USAGE OF ORTHOMCL.PL###################################
__DATA__
ORTHOMCL [2005-03-14]

Copyright (C) 2004 by University of Pennsylvania, Philadelphia, PA USA.
All rights reserved.

Before orthomcl.pl can be used, some variables (including directory variables
or parameter variables) in orthomcl_module.pm need to be set, as described in
ORTHOMCL_INSTALL.

Usage: orthomcl.pl --mode 1,2,3,4 or 5 <tagged arguments>

Modes:
~~~~~~

 1: OrthoMCL analysis from FASTA files

 2: OrthoMCL analysis based on former OrthoMCL run
    No BLAST or BLAST parsing performed.

 3: OrthoMCL analysis from user-provided BLAST result
    No BLAST performed.

 4: OrthoMCL analysis based on user-provided BPO (BLAST
    PARSE OUT) file and GG (Genome-Gene Index) file

 5: OrthoMCL analysis based on matrix of former 
    OrthoMCL run, but with LESS genomes
    
Arguments:
~~~~~~~~~~

 fa_files=<String>         Protein FASTA file names, with each file containing
                           protein sequences from one species, separated by 
                           comma. e.g. "Eco.fa,Sce.fa,Afu.fa"
 pv_cutoff=<Float>         P-Value Cutoff used in BLAST search and/or
                           identification of inparalogs and orthologs, 1e-5
                           (DEFAULT)
 pi_cutoff=<Int>           Percent Identity Cutoff <0-100> used in identification
                           of inparalogs and orthologs, 0 (DEFAULT)
 pmatch_cutoff=<Int>       Percent Match Cutoff <0-100> used in identification of
                           inparalogs and orthologs, 0 (DEFAULT)
 inflation=<Float>         Inflation value used in MCL algorithm, 1.5 (DEFAULT)
                           Increasing the inflation value increases cluster 
                           tightness, and also the number of clusters.
 former_run_dir=<String>   Former run directory, required in Mode 2, e.g. 
                           "July_21". Then the blast result file and bpo file
                           in former run directory will be used instead of 
                           running through from the very beginning.
 usr_blast_file=<String>   Blast out file provided by user, required in Mode 3
                           It will be parsed into BPO file by BioPerl
 usr_bpo_file=<String>     BPO (BLAST PARSE OUT) file provided by user, required
                           in Mode 4. About its format, please refer to 
                           ORTHOMCL_INSTALL
 usr_gg_file=<String>      GG (Genome Gene Relation) file provided by user, 
                           required in Mode 3 & 4. About its format, please refer 
                           to ORTHOMCL_INSTALL
 usr_taxa_file=<String>    TAXA file provided by user, required in Mode 5. 
                           About its format, please refer to ORTHOMCL_INSTALL
 
#### Copy and Paste from here to a BLOSUM62 file

#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 


