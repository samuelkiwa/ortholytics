# -*-Perl-*-
# launch_app.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella, Filipe Vieira
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

launch_app.PLS - Launch a system process, send an email when job finished.

=head1 SYNOPSIS

perl launch_app.PLS ./hyphymp analysis.bf inputoptfile outputfile "myemail@address.com"

=head1 DESCRIPTION

Launch a system process, send an email when job finished. Requires a
linux machine with sendmail and Mail::Sendmail. Put this script in the
same directory where binary is present.

=head1 AUTHOR - Albert Vilella

Email 

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

Filipe Vieira

=cut


# Let the code begin...

use strict;
use Mail::Sendmail;
use Getopt::Long;

my ($email,$call,$outfile,$errfile);
$errfile = "stderr";
$outfile = "stdout";

&GetOptions
    (
     "email=s"       => \$email,
     "call=s"       => \$call,
     "outfile|o=s"   => \$outfile,
     "errfile|e=s"   => \$errfile,
    );

$call =~ s/\,/\ /g;

# my (@array, @arg);
# for(my $cnt=0; $#ARGV +1 != 0; $cnt++)
# {
#     next if $arg[$cnt] =~ /\-\-/;
#     $arg[$cnt] = shift;
# }
# @array = split (m/\//, $call);
# @array = reverse(@array);

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $start_time = sprintf ("%04d%02d%02d_%02d%02d%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec);

$outfile .= ".".$start_time.".out.txt";
my $redirection = "1>$outfile";
if (-e $outfile) {unlink($outfile);}

$errfile .= ".".$start_time.".err.txt";
$redirection .= " 2>$errfile";
if (-e $errfile) {unlink($errfile);}

my $commandstring = "nohup nice $call $redirection";
printf("# $commandstring\n");

run($commandstring);

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $end_time = sprintf ("%04d%02d%02d_%02d%02d%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec);

my $message = "Command Line:\n" . $commandstring . "\n\n";
$message .= " Started at $start_time\nFinished at $end_time\n\n";
$message .= "Call:\n $call\n\n";
$message .= "Outfile:\n $outfile\n\n";
$message .= "Errfile:\n $errfile\n\n";

if(!-z $errfile)
{
    printf "Program finished but with errors: $errfile";
    $message .= "Program finished but with errors: $errfile\n\n";
}

else {unlink($errfile);}


my %mail = ( To      => "$email",
             From    => 'launch_app@localhost',
             Subject => "Sendmail batch: job finished.",
             Message => "$message"
        );

sendmail(%mail) or die $Mail::Sendmail::error;
print "OK. Log says:\n", $Mail::Sendmail::log;
print "\n";

exit(0);

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
