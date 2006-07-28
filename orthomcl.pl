#!/usr/bin/perl -w
## AUTHOR: Li Li, Feng Chen <fengchen@sas.upenn.edu>
## ORTHOMCL [2005-03-14] Version 1.2

## Copyright (C) 2004 by University of Pennsylvania, Philadelphia, PA USA.
## All rights reserved.

## Before orthomcl.pl can be used, some variables (including directory variables
## or parameter variables) in orthomcl_module.pm need to be set, as described in
## ORTHOMCL_INSTALL.

my $starttime = `date`;
use strict;
use Getopt::Long;
use File::Basename;
use orthomcl_module;

my ($mode,$fa_files,$pv_cutoff,$pi_cutoff,$pmatch_cutoff,$inflation,$maximum_weight);
my ($usr_blast_file,$usr_bpo_file,$usr_gg_file,$usr_taxa_file,$former_run_dir,$email);         # For Mode 2, 3 or 4


my $command=basename($0)." ".join(' ',@ARGV)."\n";

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
                        "email=s"             => \$email
);

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
		&constructDirectory($starttime);
		&constructAllFasta($fa_files,$all_fa_file);                                #construct all.fa file
		&executeFORMATDB($all_fa_file);                                            # and run blast
		&executeBLASTALL($all_fa_file,$blast_file,$all_fa_file,$pv_cutoff);
		&blast_parse($blast_file,$bpo_file,$pv_cutoff) unless (-e $bpo_file);
	} else {dieWithUnexpectedError("In Mode 1, NAMES OF FASTA FILES need to be given!");}
}
elsif ($mode == 2) {
	if (defined $former_run_dir) {
		&constructDirectory($starttime,$former_run_dir);
		&read_ggfile($genome_gene_file);
	} else {dieWithUnexpectedError("In Mode 2, FORMER RUN DIRECTORY needs to be given!");}
}
elsif ($mode == 3) {
	if ((defined $usr_blast_file) && (defined $usr_gg_file)) {
		&constructDirectory($starttime);
		$all_fa_file      = 'N/A';
		$genome_gene_file = $usr_gg_file;
		read_ggfile($genome_gene_file);
		$blast_file       = $usr_blast_file;
		&blast_parse($blast_file,$bpo_file,$pv_cutoff) unless (-e $bpo_file);
	} else {dieWithUnexpectedError("In Mode 3, BLAST OUT FILE and GENOME-GENE FILE are required!");}
}
elsif ($mode == 4) {
	if ((defined $usr_bpo_file) && (defined $usr_gg_file)) {
		&constructDirectory($starttime);
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
		mode5($starttime,$command,$former_run_dir,$usr_taxa_file,$inflation);
		my $endtime = `date`;
		&write_endtime_in_parameter_log($endtime);
		write_log("\nStart Time: $starttime\nEnd Time:   $endtime\n");
		die "\nStart Time: $starttime\nEnd Time:   $endtime\n";

	} else {dieWithUnexpectedError("In Mode 5, FORMER RUN DIR and TAXA LIST FILE are required!");}
} 
else {dieWithUnexpectedError("Mode 1,2,3 or 4 needs to be given!");}


&write_parameter_log($starttime,$command,$mode,$pv_cutoff,$pi_cutoff,$pmatch_cutoff,$inflation,$maximum_weight);

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

my $endtime = `date`;
print STDERR "\nStart Time: $starttime\nEnd Time:   $endtime\n";
write_log("\nStart Time: $starttime\nEnd Time:   $endtime\n");
&write_endtime_in_parameter_log($endtime);

if ($email) {
    require Mail::Sendmail;
    my $message = "Start Time: $starttime\nEnd Time:   $endtime\n";
    my %mail = ( To      => "$email",
                 From    => 'orthomcl@localhost',
                 Subject => "Sendmail batch: $0 job finished.",
                 Message => "$message"
               );

    sendmail(%mail) or die $Mail::Sendmail::error;
    print STDERR "Mail sent. Log says:\n", $Mail::Sendmail::log;
    print STDERR "\n";
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
 
