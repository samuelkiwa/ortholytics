#!/usr/local/bin/perl
# extract_feature_seq.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

extract_feature_seq - extract the corresponding sequence for a specified feature type

=head1 SYNOPSIS

extract_feature_seq.PLS -i file --format genbank --feature CDS -o output.fa [-filterproduct ribosomal]

=head1 DESCRIPTION 

This script will extract the sequence for all the features you specify.

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

 Jason Stajich <jason@bioperl.org>

=cut

use Bio::SeqIO;
use Getopt::Long;

my ($input,$format,$featuretype,$output,$filterproduct);
$featuretype ='CDS';
GetOptions(
	   'i|input:s' => \$input,
	   'format:s'  => \$format,
	   'feature:s' => \$featuretype,
	   'o|output:s'=> \$output,
           'f|filterproduct:s'=> \$filterproduct);

$input || shift if @ARGV;

my $in = new Bio::SeqIO(-file => $input,
			-format => $format);
my $out;
if ($output ) {
    $out = new Bio::SeqIO(-file => ">$output");
} else { 
    $out = new Bio::SeqIO(); # use STDOUT for output
}

my $count = 1; my $done = 0;
my $tagval = "";
print "Loading infile...\n";
while( my $seq = $in->next_seq ) {    
print "Parsing features...\n";
    foreach my $f ( grep { $_->primary_tag =~ /$featuretype/i } 
		    $seq->get_SeqFeatures ) {
	my $s = $f->spliced_seq;
        $done = 0;
	if( $featuretype =~ /gene/ ) {
	    $s->display_id($f->has_tag('gene') ? join(',',sort $f->each_tag_value('gene')) :
			   $f->has_tag('label') ? join(',',$f->each_tag_value('label')): 
                           $f->has_tag('product') ? join(',',$f->each_tag_value('product')) :
			   $s->display_id);
        } elsif ($featuretype =~ /CDS/ ) {
            if ($f->has_tag('db_xref')) {
                foreach $tagval ($f->get_tag_values('db_xref')) {
                    if(index($tagval,"GI:") == 0) {
                        unless ($filterproduct) {$s->display_id("$tagval");$done = 1}
                        last;
                    }
                }
            }
            unless ($done == 1) {
                $s->display_id($f->has_tag('gene') ? join(',',sort $f->each_tag_value('gene')) :
                               $f->has_tag('label') ? join(',',$f->each_tag_value('label')): 
                               $f->has_tag('product') ? join(',',$f->each_tag_value('product')) :
                               $s->display_id);
            }
	} else {
	    $s->display_id(sprintf("%s_%s_%d",
				   $seq->display_id, 
				   $f->primary_tag,
				   $count++));
	}

        if($filterproduct) {
            my $product = $s->display_id;
            if ($product =~ /$filterproduct/) {
	    $s->display_id(sprintf("%s_%s",
				   $tagval, 
				   $s->display_id));
                $out->write_seq($s);
            } else {
                #SKIP THIS ONE
            }
        } else {
            $out->write_seq($s);
        }
    }
}
