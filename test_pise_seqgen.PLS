use Bio::Tools::Run::AnalysisFactory::Pise;

use Getopt::Long;

my ($inputfile,$yesnoopt);

GetOptions(
	   'i|input|inputfile:s' => \$inputfile,
           'yesnoopt' => \$yesnoopt,
          );

# Build a Pise factory       
my $factory = new Bio::Tools::Run::AnalysisFactory::Pise();

# Then create an application object (Pise::Run::Tools::PiseApplication):
my $program = $factory->program('seqgen');

# Set parameters
$program->tree($inputfile);

# Next, run the program
# (notice that you can set some parameters at run time)
my $job = $program->run
    (
#     -parameter_file => "Arabidopsis.smat"
    );

# Test for submission errors:
if ($job->error) {
   print "Job submission error (",$job->jobid,"):\n";
   print $job->error_message,"\n";
   exit;
}

# Get results
print STDERR $job->content('genscan.out');
# or:
my $result_file = $job->save('genscan.out');
