# -*-Perl-*-
use strict;

# launch_psql_command.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

launch_psql_command.PLS - Launch a psql_command process, send an email when job finished.

=head1 SYNOPSIS

nice perl launch_psql_command.PLS dbname /full/path/command.sql "myemail@address.com"

=head1 DESCRIPTION

Launch a psql process, send an email when job finished. Requires a
linux machine with sendmail and Mail::Sendmail. The output of this
script process will be the same directory as specified in the
/full/path/command.sql.

=head1 AUTHOR - Albert Vilella

Email 

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut


# Let the code begin...

use strict;
use Mail::Sendmail;

my $dbname = shift;
my $command = shift;
my $email = shift;
my $outfile = "$command.stdout";
my $errfile = "$command.stderr";
my $redirection = " 1>$outfile 2>$errfile";
my $commandstring = "";
my $message;

if (-e $outfile) {unlink($outfile);}
if (-e $errfile) {unlink($errfile);}

$commandstring = "psql -d $dbname < $command $redirection";
printf("# $commandstring\n");
run($commandstring);

$message = "Output ready at $command.stdout and $command.stderr";

my %mail = ( To      => "$email",
             From    => 'launch_psql_command@localhost',
             Subject => "Sendmail batch: $0 job finished.",
             Message => "$message"
        );

unless(-z $errfile) {
    printf "Program finished but with errors: $errfile";
    $message = "Program finished but with errors: $errfile";
}

if( !-e $outfile || -z $outfile ) {
  warn( "call crashed: $? [command $commandstring]\n");
#  $message = "call crashed: $? [command $commandstring]\n";
  return undef;
}


sendmail(%mail) or die $Mail::Sendmail::error;
print "OK. Log says:\n", $Mail::Sendmail::log;

1;

################################################################################
## FUNCTIONS
################################################################################

sub run {
    my @args = @_;
    my $sysstat = system(@args);
    if ($sysstat == -1) {
        print "failed to execute: $!\n";
    }
    elsif ($sysstat & 127) {
        printf "child died with signal %d, %s coredump\n",
            ($sysstat & 127),  ($sysstat & 128) ? 'with' : 'without';
    }
    else {
        #      successfull run
        #      printf "child exited with value %d\n", $sysstat >> 8;
    }
}

1;
