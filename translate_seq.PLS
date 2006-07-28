#!perl -w
use strict;
# $Id: translate_seq.PLS,v 1.5 2005/04/26 15:08:43 jason Exp $

=head1 NAME

translate_seq - translates a sequence

=head1 SYNOPSIS

translate_seq

	   'input:s'  => \$input,
	   'output:s'  => \$output,
	   'format:s'  => \$format,


=head1 DESCRIPTION 

The script will translate one fasta file (on stdin) to protein on stdout

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

 http://bugzilla.bioperl.org/

=head1 AUTHOR

  Ewan Birney E<lt>birney@ebi.ac.ukE<gt>

=cut

use Bio::SeqIO;
use Getopt::Long;

my ($format) = 'fasta';
my ($input,$output,$tableid);

$tableid = 1;

GetOptions(
	   'format:s'  => \$format,
	   'input:s'  => \$input,
	   'output:s'  => \$output,
	   'table|tableid:s'  => \$tableid,
	   );

my $oformat = 'fasta';

# this implicity uses the <> file stream
my $seqin = Bio::SeqIO->new( -format => $format, -file => "$input"); 
my $seqout = Bio::SeqIO->new( -format => $oformat, -file => ">$output" );


while( (my $seq = $seqin->next_seq()) ) {
    my $pseq = $seq->translate(undef, undef, undef, $tableid, 'fullCDS', undef, undef);
    $seqout->write_seq($pseq);
}

__END__
