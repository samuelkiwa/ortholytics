#!/usr/local/bin/perl

# codeml_prepare_batch.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

codeml_prepare_batch.PLS - DESCRIPTION 

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Options

	   'w|work|workdir:s'  => \$workdir,
	   'e|exe|exedir:s'    => \$exedir,
	   's|seq|seqfile:s'   => \$seqfile,
	   't|tree|treefile:s' => \$treefile,
           'b|branch|branch_model|branchModel|model:s' => \$branchModel,
	   'n|ns|nssites:s'    => \$nssites,
	   'c|cod|codonfreq:s' => \$codonfreq,
	   'class|queue:s'    => \$class,
	   'mgene:s'    => \$mgene,
	   'nomail'    => \$nomail,
	   'zip|zipped'    => \$zip,
	   'sendmail' => \$sendmail,


=head1 AUTHOR - Albert Vilella

Email 

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut

# Let the code begin...
use strict;
use Getopt::Long;
use Bio::Tools::Run::Phylo::PAML::Codeml;
use File::Path;

my ($workdir, $exedir, $seqfile, $treefile, $branchModel, $nssites,
    $codonfreq, $class, $mgene, $sendmail, $nomail, $zip, $verbose,
    $startmail, $rateancestor, $ncatg) =
    undef;

my ($runid,$program) = undef;
$branchModel = 0; $nssites = 0; $codonfreq = 5;
$runid = 1;
$program = "codeml";
$class = "research_long";
$mgene = 0;
$verbose = 0;
$rateancestor = 0;
$ncatg = 3;

GetOptions(
	   'w|work|workdir:s'  => \$workdir,
	   'e|exe|exedir:s'    => \$exedir,
	   's|seq|seqfile:s'   => \$seqfile,
	   't|tree|treefile:s' => \$treefile,
           'b|branch|branch_model|branchModel|model:s' => \$branchModel,
	   'n|ns|nssites:s'    => \$nssites,
	   'c|cod|codonfreq:s' => \$codonfreq,
	   'class|queue:s'    => \$class,
	   'mgene:s'    => \$mgene,
	   'ncatg:s'    => \$ncatg,
	   'zip|zipped'    => \$zip,
           'verbose:s' => \$verbose,
           'a|rateancestor:s' => \$rateancestor,
	   'nomail'    => \$nomail,
	   'sendmail' => \$sendmail,
	   'startmail'    => \$startmail,
          );

my $analysis;
my ($dir, $batch, $loader);
$runid = sprintf("%04d", $runid);
$analysis = "$program"."$runid";
$dir = "$workdir/$analysis";
while (-d $dir) {
    $runid++;
    $runid = sprintf("%04d", $runid);
    $analysis = "$program"."$runid";
    $dir = "$workdir/$analysis";
}

unless (-d $dir) {
    File::Path::mkpath($dir);
#     mkdir $dir or die "couldnt create directory: $!";
}


my $alignio = new Bio::AlignIO
    (
     -format => "phylip",
     -file   => "$seqfile",
    );

my $treeio = new Bio::TreeIO
    (
     -format => 'newick', 
     -file => "$treefile",
    );

# my $tree = $treeio->next_tree;
# my $aln = $alignio->next_aln;

my $codeml = new Bio::Tools::Run::Phylo::PAML::Codeml();
# $codeml->alignment($aln);
# $codeml->tree($tree);

$codeml->no_param_checks(1);
if ($branchModel =~ /m0/i) {$branchModel = "0";} 
elsif ($branchModel =~ /fr/i) {$branchModel = "1";}
$codeml->set_parameter("model","$branchModel");
$codeml->set_parameter("noisy","9");
$codeml->set_parameter("verbose","$verbose");
$codeml->set_parameter("runmode","0");
$codeml->set_parameter("seqtype","1");
$codeml->set_parameter("CodonFreq","$codonfreq");
$codeml->set_parameter("clock","0");
$codeml->set_parameter("aaDist","0");
$codeml->set_parameter("aaRatefile","jones.dat");
$codeml->set_parameter("NSsites","$nssites");
$codeml->set_parameter("icode","0");
$codeml->set_parameter("Mgene","$mgene");
$codeml->set_parameter("fix_kappa","0");
$codeml->set_parameter("kappa","2");
$codeml->set_parameter("fix_omega","0");
$codeml->set_parameter("fix_alpha","1");
$codeml->set_parameter("alpha","0.");
$codeml->set_parameter("Malpha","0");
$codeml->set_parameter("ncatG","$ncatg");
$codeml->set_parameter("getSE","0");
$codeml->set_parameter("RateAncestor","$rateancestor");
$codeml->set_parameter("Small_Diff",".5e-6");
$codeml->set_parameter("cleandata","1");
$codeml->set_parameter("fix_blength","1");
$codeml->set_parameter("method","0");
$codeml->save_tempfiles(1);

$codeml->tempdir($dir);
$codeml->prepare($seqfile,$treefile);
print "# Creating $dir\n";
print " nice -n 19 sh $dir/$analysis.loader &\n";
print "# bsub < $dir/$analysis.loader\n";

$loader = "$workdir/$analysis/$analysis.loader";
open (LOADERFILE, ">$loader") or die "cannot open file $loader: $!";

print LOADERFILE "# BSUB -J $analysis\n";
print LOADERFILE "# BSUB -q $class\n";
print LOADERFILE "# BSUB -o $workdir/$analysis/stdout.$analysis.log\n";
print LOADERFILE "# BSUB -e $workdir/$analysis/stderr.$analysis.log\n";
print LOADERFILE "# BSUB -B -u avilella\@gmail.com\n" if(!$nomail && $startmail);
print LOADERFILE "# BSUB -N -u avilella\@gmail.com\n" unless($nomail);
print LOADERFILE "\n";
print LOADERFILE "cd $dir\n";
print LOADERFILE "echo\n";
print LOADERFILE "\ndate +2%3y%m%d_%H%M%S\n\n";
print LOADERFILE "\necho \`date +2%3y%m%d_%H%M%S\` > $workdir/$analysis/xpent.txt\n\n";
print LOADERFILE "$exedir/$program codeml.ctl > $dir/stdout.log \n";
print LOADERFILE "rm rst* \n" unless ($verbose);
print LOADERFILE "rm 2N* \n" unless ($verbose>1);
print LOADERFILE "rm lnf \n" unless ($verbose>1);
print LOADERFILE "rm rub \n" unless ($verbose>1);
print LOADERFILE "bzip2 --best mlc \n" if ($zip);
print LOADERFILE "bzip2 --best *.log \n" if ($zip);
print LOADERFILE "\necho \`date +2%3y%m%d_%H%M%S\` >> $workdir/$analysis/xpent.txt\n\n";
print LOADERFILE "\ndate +2%3y%m%d_%H%M%S\n\n";
print LOADERFILE "\nperl ~/ortholytics/launch_app.pl --call \"echo\,\'$workdir/$analysis/$analysis.loader\'\" -email \"avilella\@gmail.com\"\n" if ($sendmail);

close LOADERFILE;

1;
