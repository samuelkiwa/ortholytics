#!/usr/local/bin/perl
# phylo_prepare_batch.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

phylo_prepare_batch.PLS - DESCRIPTION 

=head1 SYNOPSIS

phylo_prepare_batch.PLS 
-workdir /my/work/dir         # the dir where the work run will be saved
-exedir /my/exe/dir           # the dir where the exe file is located
-seqfile file.phy
-treefile file.tree

=head1 DESCRIPTION

This script will create a set of directories for a set of phylogenetic
analyses one wants to conduct.

It will create, but not execute, the dirs with the ready-to-run
elements of the executables according to the various input data files
and parameters.

type_ana - type of analysis to run:

-a 1 or /AnalyseCodonData/i [default]

-a 2 /NucModelCompare/i

-a 3 /CodonModelCompare/i

-a 4 /Bivariate/i

-a 5 /FitIndividualDNs/i

-a 6 /gabranch/i

Options:

	   'w|work|workdir:s' => \$workdir,
	   'tmp|tempdir:s' => \$tmpdir,
	   'e|exe|exedir:s' => \$exedir,
	   's|seq|seqfile:s' => \$seqfile,
	   'g|gl|genelenfile:s' => \$genelenfile,
	   't|tree|treefile:s' => \$treefile,
	   'c|cpu|cpunum:s' => \$cpunum,
	   'class:s' => \$class,
	   'i|iter|iterations:s' => \$iterations,
	   'pboot' => \$pboot,
           'a|type|type_ana|ana|analysis:s' => \$type_ana,
           'b|branch|branch_model|branchModel|model:s' => \$branchModel,
           'cf|codon|codonfreq:s' => \$codFreqModel,
           'd|dist|distrib|distribution:s' => \$distrib,
           'r|rate|rateclas:s' => \$rateclasses,
	   'ncatg:s'    => \$ncatg,
           'dndspg:s' = \$dndspg;
           'bfpg|basefreqpg:s' = \$basefreqpg;
           'nrbpg|nuclratebiaspg:s' = \$nuclratebiaspg;
           'blpg|branchlengthpg:s' = \$branchlengthpg;
           'opt|optprec|prec:s' => \$optprec,
           'l|load|loader:s' => \$type_loader,
           'q|queue|queuename:s' => \$queue,
           'p|program|binary:s' => \$program,
	   'zip|zipped' => \$zip,
           'sendmail' => \$sendmail,
           'nomail' => \$nomail,
           'startmail' => \$startmail,


dNdS pg
"Shared","All genes share the same dN/dS.",
"By gene","Each gene has its own dN/dS.");

Base Frequencies pg
"Shared","All genes share the same base frequencies.",
"By gene","Each gene has its own base frequencies");

Nucleotide Biases pg
"Shared","All genes share the same nucleotide rate bias parameters.",
"By gene","Each gene has its own nucleotide rate bias parameters");

Branch Length pg
"Shared","All genes share the same branch length (may not be implemented for all options).",
"By gene","Each gene has its own branch lengths",
"Proportional","Branch lengths are proportional between genes");


=head1 AUTHOR - Albert Vilella

Email 

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut


# Let the code begin...
use strict;
use Getopt::Long;
use File::Path;

my $default_email = "avilella\@gmail.com"; # PLEASE change!!!

my ($workdir, $tmpdir, $exedir, $seqfile, $treefile, $genelenfile,
    $cpunum, $class, $iterations, $pboot, $type_ana, $branchModel,
    $codFreqModel, $distrib, $rateclasses, $type_loader, $queue,
    $dndspg, $basefreqpg, $nuclratebiaspg, $branchlengthpg, $noconstraint,
    $optprec, $program, $sendmail, $nomail, $startmail, $zip, $ncatg, 
    $gasitevar, $gabclasses, $garates);

$dndspg = 2;
$basefreqpg = 1;
$nuclratebiaspg = 1;
$branchlengthpg = 2;
$gasitevar = 2;
$gabclasses = 2;
$garates = 3;

my ($runid) = undef;
$cpunum = 1;
$runid = 1;
$class = "medium";
$type_ana = 1;
$iterations = 0;
$branchModel = 1;
$codFreqModel = "MG94CUSTOM";
$type_loader = 1;
$program = "HYPHYMP";
$queue = "research_long";
$distrib = 1;
$rateclasses = 3;
$optprec="0.00001";
#$maximum_iterations_per_variable="999999999";
$ncatg = 3;

GetOptions(
	   'w|work|workdir:s' => \$workdir,
	   'tmp|tempdir:s' => \$tmpdir,
	   'e|exe|exedir:s' => \$exedir,
	   's|seq|seqfile:s' => \$seqfile,
	   'g|gl|genelenfile:s' => \$genelenfile,
	   't|tree|treefile:s' => \$treefile,
	   'c|cpu|cpunum:s' => \$cpunum,
	   'class:s' => \$class,
	   'i|iter|iterations:s' => \$iterations,
	   'pboot' => \$pboot,
           'a|type|type_ana|ana|analysis:s' => \$type_ana,
           'b|branch|branch_model|branchModel|model:s' => \$branchModel,
           'cf|codon|codonfreq:s' => \$codFreqModel,
           'd|dist|distrib|distribution:s' => \$distrib,
           'r|rate|rateclas:s' => \$rateclasses,
	   'ncatg:s'    => \$ncatg,
           'dndspg:s' => \$dndspg,
           'bfpg|basefreqpg:s' => \$basefreqpg,
           'nrbpg|nuclratebiaspg:s' => \$nuclratebiaspg,
           'blpg|branchlengthpg:s' => \$branchlengthpg,
           'noconstraint' => \$noconstraint,
           'gasitevar:s' => \$gasitevar,
           'gabclasses:s' => \$gabclasses,
           'garates:s' => \$garates,
           'opt|optprec|prec:s' => \$optprec,
           'l|load|loader:s' => \$type_loader,
           'q|queue|queuename:s' => \$queue,
           'p|program|binary:s' => \$program,
	   'zip|zipped' => \$zip,
           'sendmail' => \$sendmail,
           'nomail' => \$nomail,
           'startmail' => \$startmail,
          );

my $nodes = $cpunum/2;
if ($branchModel =~ /m0/i) {$branchModel = "2";} 
elsif ($branchModel =~ /fr/i) {$branchModel = "1";}

# Choose which set
my ($AnalyseCodonData, $NucModelCompare, 
    $CodonModelCompare, $Bivariate, $FitIndividualDNs);
# common vars
my ($chooseAna, $geneticCode, $modelRejection, 
    $summaryResults, $quitResults, $startResults, 
    $resultsFile, $post_ana);

########################################

# Directory creation
my $analysis;
my ($dir, $batch, $loaderfile);

$runid = sprintf("%04d", $runid);
my $hyphybinname = lc("$program");
$analysis = "$hyphybinname"."$runid";
$dir = "$workdir/$analysis";
while (-d $dir) {
    $runid++;
    $runid = sprintf("%04d", $runid);
    $analysis = "$hyphybinname"."$runid";
    $dir = "$workdir/$analysis";
}
unless (-d $dir) {
    File::Path::mkpath($dir);
}

########################################
# Analysis Type

# AnalyseCodonData

my $noboot = 0;
$noboot = 1 if ($iterations == 0);
$chooseAna = "1\n1\n";
$geneticCode = "1\n";
$branchModel .= "\n";
$codFreqModel .= "\n";
my $nuclModel = "";
if ($codFreqModel =~ /MG94CUSTOM/) {
    $nuclModel = "010010\n";
}
if ($branchModel =~ /[3-4]/) {
    $distrib .= "\n";
    $branchModel .= $distrib;
    $rateclasses .= "\n";
    $branchModel .= $rateclasses;
    $treefile .= "\n2";
    $noboot = 1;
}
$startResults = "";
$startResults .= "y\n" unless ($branchModel =~ /[3-4]/);

my $resoptnum;
$resoptnum = 2;
$summaryResults = 
    "y\n$resoptnum\n1\n" . "$dir/res.1.txt\n". 
    "y\n$resoptnum\n2\n" . "$dir/res.2.txt\n". 
    "y\n$resoptnum\n3\n" . "$dir/res.3.txt\n". 
    "y\n$resoptnum\n4\n" . "$dir/res.4.txt\n". 
    "y\n$resoptnum\n5\n" . "$dir/res.5.txt\n". 
    "y\n$resoptnum\n6\n" . "$dir/res.6.txt\n";
$summaryResults .= "y\n3\n4\n" if ($branchModel =~ /[3-4]/);
$resoptnum = 8; $resoptnum++ if ($branchModel =~ /[3-4]/);
$summaryResults .= 
    "y\n$resoptnum\n2\n" . "$dir/res.dSdNtrees.txt\n";

#Variance Estimates
my $varest;
$resoptnum = 3;
$varest .= 
    "y\n$resoptnum\n3\ny\n$resoptnum\n2\ny\n$resoptnum\n1\n";
$resoptnum = 8;
$varest .= 
    "y\n$resoptnum\n1\n"; 

$resoptnum = 6;
my $Vest_NParBoot = 
    "y\n$resoptnum\n" . 
    "$iterations\n" . 
    "$dir/res.Vest." . 
    "$iterations" . 
    ".npboot.csv\n" . 
    "$varest";

$resoptnum = 5; 
my $Vest_ParBoot = 
    "y\n$resoptnum\n" . 
    "$iterations\n" . 
    "$dir/res.Vest." . 
    "$iterations" . 
    ".pboot.csv\n" . 
    "$varest";

$quitResults = "n\n";

$AnalyseCodonData = 
    "$chooseAna" . 
    "$geneticCode" . 
    "$seqfile\n" . 
    "$codFreqModel" . 
    "$branchModel" . 
    "$nuclModel" . 
    "$treefile\n" . 
    "$startResults" . 
    "$summaryResults";
$AnalyseCodonData .=
    "$Vest_ParBoot" unless ($noboot);
$AnalyseCodonData .=
    "$Vest_NParBoot" if ((!$noboot) && (!$pboot));
$AnalyseCodonData .=
    "$quitResults";

my $postAnaAnalyseCodonData = undef;

########################################

# NucModelCompare

$chooseAna = "6\n4\n";
$branchModel = "1\n";
$modelRejection = "0.05\n";
$resultsFile = "$dir/res.NucModelCompare.txt\n";

$summaryResults = 
    "y\n2\n1\n" . "$dir/res.1.txt\n". 
    "y\n2\n2\n" . "$dir/res.2.txt\n". 
    "y\n2\n3\n" . "$dir/res.3.txt\n". 
    "y\n2\n4\n" . "$dir/res.4.txt\n". 
    "y\n2\n5\n" . "$dir/res.5.txt\n". 
    "y\n2\n6\n" . "$dir/res.6.txt\n";

$quitResults = "n\n";

$NucModelCompare = 
    "$chooseAna" . 
    "$branchModel" . 
    "$seqfile\n" . 
    "$treefile\n" . 
    "$modelRejection" . 
    "$resultsFile" . 
    "$summaryResults" . 
    "$quitResults";

my $deleteModelFiles = 
    "rm -f $dir/res.NucModelCompare.txt.*\n";

my $postAnaNucModelCompare =
    "$deleteModelFiles";

########################################

# CodonModelCompare

$chooseAna = "6\n1\n";
$geneticCode = "1\n";
$branchModel .= "\n";
my $estimateBranchLengths = "1\n";
$modelRejection = "0.05\n";
$resultsFile = "$dir/res.CodonModelCompare.txt\n";

$summaryResults = 
    "y\n2\n1\n" . "$dir/res.1.txt\n". 
    "y\n2\n2\n" . "$dir/res.2.txt\n". 
    "y\n2\n3\n" . "$dir/res.3.txt\n". 
    "y\n2\n4\n" . "$dir/res.4.txt\n". 
    "y\n2\n5\n" . "$dir/res.5.txt\n". 
    "y\n2\n6\n" . "$dir/res.6.txt\n";

$quitResults = "n\n";

$CodonModelCompare = 
    "$chooseAna" . 
    "$geneticCode" . 
    "$seqfile\n" . 
    "$treefile\n" . 
    "$branchModel" . 
    "$estimateBranchLengths" . 
    "$modelRejection" . 
    "$resultsFile" . 
    "$summaryResults" . 
    "$quitResults";

my $postAnaCodonModelCompare = undef;

########################################

# Bivariate

$chooseAna = "2\n2\n1\n";
$geneticCode = "1\n";
$branchModel .= "\n";
$codFreqModel .= "\n";
my $nuclModel = "";
if ($codFreqModel =~ /MG94CUSTOM/) {
    $nuclModel = "010010\n";
}
my $joinly_opt = "1\n";
my $initial_val = "1\n";
my $strat = $ncatg-2;
my $stratification;
if (!$noconstraint) {
    $stratification = "2\n$strat\n1\n";
} else {
    $stratification = "1\n";
}
my ($sec,$min,$hour,$mday,$mon,
    $year,$wday,$yday,$isdst) = localtime(time);
my $timestamp = sprintf
    ("%04d%02d%02d_%02d%02d%02d", 
     $year+1900,$mon+1,$mday,$hour,$min,$sec
    );
my $rand = int(rand(99999));
my $resfile = "res.bivariate.$ncatg.$timestamp.$rand";
my $res_nuc_lf = "y\n1\n1\n6\n";
my $res_lf = "y\n1\n2\n6\n";
$quitResults = "n\n";

$Bivariate = 
    "$chooseAna" . 
    "$seqfile\n" . 
    "$geneticCode" . 
    "$treefile\n" . 
    "$joinly_opt" . 
    "$nuclModel" . 
    "$initial_val" . 
    "$ncatg\n" . 
    "$stratification" . 
    "$tmpdir/$resfile\n" . 
    "$res_nuc_lf" . 
    "$res_lf";
$Bivariate .=
    "$quitResults";

my $deleteModelFiles = "rm $tmpdir/$resfile\n";
my $moveresfiles .= 
    "mv $tmpdir/res.bivariate.$ncatg.$timestamp.$rand.summary $dir/res.bivariate.$ncatg.summary.txt\n" . 
    "mv $tmpdir/res.bivariate.$ncatg.$timestamp.$rand.posterior $dir/res.bivariate.$ncatg.posterior.csv\n" . 
    "mv $tmpdir/res.bivariate.$ncatg.$timestamp.$rand.cAIC $dir/res.bivariate.$ncatg.cAIC.txt\n";

my $postAnaBivariate =
    "$deleteModelFiles" . 
    "$moveresfiles";

########################################

# FitIndividualDNs

$geneticCode = "1\n";
$codFreqModel .= "\n";
my $nuclModel = "";
if ($codFreqModel =~ /MG94CUSTOM/) {
    $nuclModel = "010010\n";
}
my $dndspergene = $dndspg;
my $basefreqpergene = $basefreqpg;
my $nuclratebiaspergene = $nuclratebiaspg;
my $branchlengthpergene = $branchlengthpg;
my $outputoption = "1\n";
$resoptnum = 2;
$summaryResults = 
    "y\n$resoptnum\n1\n" . "$dir/res.1.txt\n". 
    "y\n$resoptnum\n2\n" . "$dir/res.2.txt\n". 
    "y\n$resoptnum\n3\n" . "$dir/res.3.txt\n". 
    "y\n$resoptnum\n4\n" . "$dir/res.4.txt\n". 
    "y\n$resoptnum\n5\n" . "$dir/res.5.txt\n". 
    "y\n$resoptnum\n6\n" . "$dir/res.6.txt\n";
$varest = "";
$varest .= 
    "y\n3\n3\ny\n3\n2\ny\n3\n1\n";
my $paramsampler = 
    "y\n5\n1\n" . 
    "$dir/res.paramsampler.2.csv\n" . 
    "2\n$iterations\n$iterations\n";
my $ancestralseqs = 
    "y\n7\n" . 
    "$dir/res.ancestralseqs.fasta.log\n";
$quitResults = "n\n";

$FitIndividualDNs = 
    "$geneticCode" . 
    "$seqfile\n" . 
    "$genelenfile\n" . 
    "$nuclModel" . 
    "$treefile\n" . 
    "$dndspergene\n" . 
    "$basefreqpergene\n" . 
    "$nuclratebiaspergene\n" . 
    "$branchlengthpergene\n" . 
    "$outputoption" . 
    "$summaryResults" . 
    "$varest";
$FitIndividualDNs .= "$paramsampler" if (0 != $iterations);
$FitIndividualDNs .= 
    "$ancestralseqs";
$FitIndividualDNs .= 
    "$quitResults";

my $postAnaFitIndividualDNs = "";

########################################

# gabranch

$geneticCode = "1\n";
$codFreqModel .= "\n";
my $nuclModel = "";
if ($codFreqModel =~ /MG94CUSTOM/) {
    $nuclModel = "010010\n";
}

my $gasitevaropt = $gasitevar;
my $gabclassesopt = $gabclasses;
my $gadistribopt = "8";
my $garatesopt = $garates;
my $resfile = "$dir/res.ga.s$gasitevar.c$gabclasses.r$garates.txt";

my $gabranch = 
    "$gabclassesopt\n" . 
    "$geneticCode" . 
    "$seqfile\n" . 
    "$gasitevaropt\n" . 
    "$gadistribopt\n" . 
    "$garatesopt\n" . 
    "$nuclModel" . 
    "$resfile\n" . 
    "$treefile\n" . 
    "$quitResults";

my $postAnagabranch = "";

########################################

my $ana;
if (1 == $type_ana || $type_ana =~ /analy[sz]ecodondata/i) {
    $ana = $AnalyseCodonData;
    $post_ana = $postAnaAnalyseCodonData;
} elsif (2 == $type_ana || $type_ana =~ /nuc\S*modelcompare/i) {
    $ana = $NucModelCompare;
    $post_ana = $postAnaNucModelCompare;
} elsif (3 == $type_ana || $type_ana =~ /codonmodelcompare/i) {
    $ana = $CodonModelCompare;
    $post_ana = $postAnaCodonModelCompare;
} elsif (4 == $type_ana || $type_ana =~ /bivariate/i) {
    $ana = $Bivariate;
    $post_ana = $postAnaBivariate;
} elsif (5 == $type_ana || $type_ana =~ /FitIndividualDNs/i) {
    $ana = $FitIndividualDNs;
    $post_ana = $postAnaFitIndividualDNs;
} elsif (6 == $type_ana || $type_ana =~ /gabranch/i) {
    $ana = $gabranch;
    $post_ana = $postAnagabranch;
}

########################################
# Type of loader

print "Creating $workdir/$analysis\n\n";
print "# nice -n 19 sh $workdir/$analysis/*.loader &\n",
      "# bsub < $workdir/$analysis/*.loader \n" if (1 == $type_loader || $type_loader =~ /seq/i);
print "# llsubmit $workdir/$analysis/*.loader \n" if (2 == $type_loader || $type_loader =~ /level/i || 4 == $type_loader);
print "# bsub < $workdir/$analysis/*.loader \n" if (3 == $type_loader || $type_loader =~ /lsf/i);

$batch = "$workdir/$analysis/$analysis.bf";
open (OUTFILE,">$batch") or die "cannot open file $batch: $!";
print OUTFILE "$ana";
close OUTFILE;

$loaderfile = "$workdir/$analysis/$analysis.loader";
open (LOADERFILE, ">$loaderfile") or die "cannot open file $loaderfile: $!";

my ($sequential_run, $sequential_call, 
    $loadleveler_mpi_run, $loadleveler_mpi_call, 
    $lsf_mpi_run, $lsf_mpi_call, 
    $loadleveler_seq_run, $loadleveler_seq_call) = "";

$sequential_run .= "#### job scheduler\n";
$sequential_run .= "# BSUB -J $analysis\n";
$sequential_run .= "# BSUB -q $queue\n";
$sequential_run .= "# BSUB -o $workdir/$analysis/stdout.$analysis.log\n";
$sequential_run .= "# BSUB -e $workdir/$analysis/stderr.$analysis.log\n";
$sequential_run .= "# BSUB -B -u $default_email\n" if(!$nomail && $startmail);
$sequential_run .= "# BSUB -N -u $default_email\n" unless($nomail);
$sequential_run .= "#### job scheduler end\n\n";

$sequential_call .= 
    "\n## hyphy options\n".
    "export OPTIMIZATION_PRECISION=$optprec\n" . 
    "export BASEPATH=$exedir\n";
$sequential_call .= "export MESSAGE_LOGGING=0.0\n" if (5 == $type_ana);
$sequential_call .= "## hyphy options end\n\n";
$sequential_call .= "cd $exedir\n";
$sequential_call .= "echo\n";
$sequential_call .= "\necho \`date +2%3y%m%d_%H%M%S\` > $workdir/$analysis/xpent.txt\n\n";
$sequential_call .= "$exedir/$program CPU=$cpunum -p ";
$sequential_call .= "$exedir/MultiGene/BatchFiles/FitIndividualDNs.bf " if (5 == $type_ana);
$sequential_call .= "$exedir/gabranch/BatchFiles/BranchSelector.bf " if (6 == $type_ana);
$sequential_call .= "< $batch > $batch.stdout.log\n";
$sequential_call .= "cp $exedir/messages.log $workdir/$analysis/\n" if (5!= $type_ana);
$sequential_call .= 
    "\nperl ~/ortholytics/launch_app.pl --call \"echo\,\'$workdir/$analysis/$analysis.loader\'\" -email \"$default_email\"\n" if ($sendmail);
$sequential_call .= "\necho \`date +2%3y%m%d_%H%M%S\` >> $workdir/$analysis/xpent.txt\n\n";

##########

$loadleveler_mpi_run .= "#! /bin/tcsh\n#\n";
$loadleveler_mpi_run .= "\n#### job scheduler\n";
$loadleveler_mpi_run .= "# @ job_type = parallel\n";
$loadleveler_mpi_run .= "# @ class = $class\n";
$loadleveler_mpi_run .= "# @ output = $workdir/$analysis/hyphy.\$(jobid).out\n";
$loadleveler_mpi_run .= "# @ error =  $workdir/$analysis/hyphy.\$(jobid).err\n";
$loadleveler_mpi_run .= "# @ restart = no\n";
$loadleveler_mpi_run .= "# @ requirements = (Feature == \"myrinet\")\n";
$loadleveler_mpi_run .= "# @ total_tasks = $cpunum\n";
$loadleveler_mpi_run .= "# @ node = $nodes\n";
$loadleveler_mpi_run .= "# @ group = ub78\n";
$loadleveler_mpi_run .= "# @ environment = AUTO_PARALLELIZE_OPTIMIZE=1;\n";
#$loadleveler_mpi_run .= "# @ environment = MPIOPTIMIZER=1;\n";
$loadleveler_mpi_run .= "# @ queue\n\n";
$loadleveler_mpi_run .= "#environment\n";
$loadleveler_mpi_run .= "setenv MP_EUILIB gm\n";
$loadleveler_mpi_run .= "setenv OBJECT_MODE 64\n";
$loadleveler_mpi_run .= "setenv MP_RSH ssh\n";
$loadleveler_mpi_run .= "#Conseguimos la Machine_list del loadleveler\n";
$loadleveler_mpi_run .= "cd $exedir\n";
$loadleveler_mpi_run .= "setenv MLIST machine_list_\$LOADL_STEP_ID\n";
$loadleveler_mpi_run .= "/opt/ibmll/LoadL/full/bin/ll_get_machine_list > \$MLIST\n";
$loadleveler_mpi_run .= "set NPROCS = \`cat \$MLIST \|wc -l\`\n";
$loadleveler_mpi_run .= "\n#### job scheduler end\n\n";

#$loadleveler_mpi_call .= "setenv MPIOPTIMIZER 1\n";
$loadleveler_mpi_call .= 
    "\n## hyphy options\n".
    "setenv OPTIMIZATION_PRECISION $optprec\n" . 
    "setenv AUTO_PARALLELIZE_OPTIMIZE 1\n";
$loadleveler_mpi_call .= "setenv MESSAGE_LOGGING 0.0\n" if (5 == $type_ana);
$loadleveler_mpi_call .= "## hyphy options end\n\n";
$loadleveler_mpi_call .= "setenv MPIPARTITIONS $cpunum\n" if (5 == $type_ana);
$loadleveler_mpi_call .= "\necho \`date +2%3y%m%d_%H%M%S\` > $workdir/$analysis/xpent.txt\n\n";
$loadleveler_mpi_call .= "time mpirun -s -np \${NPROCS} -machinefile \$MLIST $exedir/$program -p ";
$loadleveler_mpi_call .= "$exedir/MultiGene/BatchFiles/FitIndividualDNs.bf " if (5 == $type_ana);
$loadleveler_mpi_call .= "$exedir/gabranch/BatchFiles/BranchSelector.bf " if (6 == $type_ana);
$loadleveler_mpi_call .= "< $batch > $batch.stdout.log\n";
$loadleveler_mpi_call .= "cp $exedir/messages.log $workdir/$analysis/\n" if (5!= $type_ana);
$loadleveler_mpi_call .= "echo \"End of execution with \$NPROCS processors\"\n";
$loadleveler_mpi_call .= "\necho \`date +2%3y%m%d_%H%M%S\` >> $workdir/$analysis/xpent.txt\n\n";

##########

$loadleveler_seq_run .= "#! /bin/tcsh\n#\n";
$loadleveler_seq_run .= "\n#### job scheduler\n";
$loadleveler_seq_run .= "#@ initialdir = $dir\n";
$loadleveler_seq_run .= "#@ job_type = serial\n";
$loadleveler_seq_run .= "#@ job_name = $analysis\n";
$loadleveler_seq_run .= "#@ output = $workdir/$analysis/stdout.$analysis.log\n";
$loadleveler_seq_run .= "#@ error = $workdir/$analysis/stderr.$analysis.log\n";
$loadleveler_seq_run .= "#@ environment = MP_SHARED_MEMORY=yes;MP_INFOLEVEL=3;MP_EUILIB=ip;MP_EUIDEVICE=en0.0\n" if (1 < $cpunum);
$loadleveler_seq_run .= "#@ class = $class\n";
$loadleveler_seq_run .= "#@ queue\n";
$loadleveler_seq_run .= "\n#### job scheduler end\n\n";
$loadleveler_seq_run .= "cd $exedir\n";

$loadleveler_seq_call .= 
    "\n## hyphy options\n" . 
    "setenv OPTIMIZATION_PRECISION $optprec\n" . 
    "setenv AUTO_PARALLELIZE_OPTIMIZE 1\n";
$loadleveler_seq_call .= "setenv MESSAGE_LOGGING 0.0\n" if (5 == $type_ana);
$loadleveler_seq_call .= "## hyphy options end\n\n";
$loadleveler_seq_call .= "\necho \`date +2%3y%m%d_%H%M%S\` > $workdir/$analysis/xpent.txt\n\n";
$loadleveler_seq_call .= "$exedir/$program -p ";
$loadleveler_seq_call .= "$exedir/MultiGene/BatchFiles/FitIndividualDNs.bf " if (5 == $type_ana);
$loadleveler_seq_call .= "$exedir/gabranch/BatchFiles/BranchSelector.bf " if (6 == $type_ana);
$loadleveler_seq_call .= "< $batch > $batch.stdout.log\n";
$loadleveler_seq_call .= "cp $exedir/messages.log $workdir/$analysis/\n" if (5!= $type_ana);
$loadleveler_seq_call .= "\necho \`date +2%3y%m%d_%H%M%S\` >> $workdir/$analysis/xpent.txt\n\n";


##########

$lsf_mpi_run .= "\n#### job scheduler\n";
$lsf_mpi_run .= "# LSBATCH: User input\n";
$lsf_mpi_run .= "#!/bin/ksh\n";
#$lsf_mpi_run .= "export MPIOPTIMIZER=1\n";
$lsf_mpi_run .= 
    "\n## hyphy options\n".
    "export OPTIMIZATION_PRECISION=$optprec\n".
    "export AUTO_PARALLELIZE_OPTIMIZE=1\n";
$lsf_mpi_run .= "export MESSAGE_LOGGING=0.0\n" if (5 == $type_ana);
$lsf_mpi_run .= "## hyphy options end\n\n";
$lsf_mpi_run .= "export MPIPARTITIONS=$cpunum\n" if (5 == $type_ana);
$lsf_mpi_run .= "# BSUB -J $analysis\n";
$lsf_mpi_run .= "# BSUB -q parallel\n";
$lsf_mpi_run .= "# BSUB -n $cpunum\n";
$lsf_mpi_run .= "# BSUB -R \"obacs\"\n";
$lsf_mpi_run .= "# BSUB -o $workdir/$analysis/stdout.$analysis.log\n";
$lsf_mpi_run .= "# BSUB -e $workdir/$analysis/stderr.$analysis.log\n";
$lsf_mpi_run .= "# BSUB -B -u $default_email\n" unless($nomail);
$lsf_mpi_run .= "# BSUB -N -u $default_email\n" unless($nomail);
$lsf_mpi_run .= "#### job scheduler end\n\n";

$lsf_mpi_call .= "\ndate +2%3y%m%d_%H%M%S\n\n";
$lsf_mpi_call .= "cd $exedir\n";
$lsf_mpi_call .= "mpirun -np $cpunum $exedir/$program -p ";
$lsf_mpi_call .= "$exedir/MultiGene/BatchFiles/FitIndividualDNs.bf " if (5 == $type_ana);
$lsf_mpi_call .= "$exedir/gabranch/BatchFiles/BranchSelector.bf " if (6 == $type_ana);
$lsf_mpi_call .= "< $batch > $batch.stdout.log\n";
$lsf_mpi_call .= "cp $exedir/messages.log $workdir/$analysis/\n" if (5!= $type_ana);
$lsf_mpi_call .= "\ndate +2%3y%m%d_%H%M%S\n\n";

########################################

# Write bf and loader file

if (1 == $type_loader || $type_loader =~ /seq/i) {
    $type_loader = "$sequential_run";
    $type_loader .= "$sequential_call";
} elsif (2 == $type_loader || $type_loader =~ /level/i) {
    $type_loader = "$loadleveler_mpi_run";
    $type_loader .= "$loadleveler_mpi_call";
} elsif (3 == $type_loader || $type_loader =~ /lsf/i) {
    $type_loader = "$lsf_mpi_run";
    $type_loader .= "$lsf_mpi_call";
} elsif (4 == $type_loader) {
    $type_loader = "$loadleveler_seq_run";
    $type_loader .= "$loadleveler_seq_call";
}

print LOADERFILE "$type_loader";
print LOADERFILE "$post_ana";
print LOADERFILE "\nbzip2 --best $dir/*.log\n" if ($zip);
print LOADERFILE "\nbzip2 --best $dir/*.csv\n" if ($zip);
close LOADERFILE;

1;
