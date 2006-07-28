#! /usr/bin/perl -w

=head1 NAME 

obtain_genbank.pl

=head1 SYNOPSIS

   obtain_genbank.pl --type acc --outformat genbank AE002122

=head1 DESCRIPTION

This script allows the dynamic retrieval of sequences from the GenBank
database at NCBI, via an Entrez query.  The files are assumed
formatted identically with the format given in the --outformat flag.

WARNING: Please do NOT spam the Entrez web server with multiple
requests.  NCBI offers Batch Entrez for this purpose.

=head1 ARGUMENTS

When querying for GenBank identifiers like 6899192 will need to call:

   obtain_genbank.pl --type id --outformat genbank 6899192

When querying for GenBank accessions like AE002122 will need to call:

   obtain_genbank.pl --type acc --outformat genbank AE002122

Default values for each parameter are shown in square brackets.

=over 2 

=item --type [id]

This if for retrieving via either identifiers "id" or accessions
"acc".

=item --outformat [genbank]

Examples: 
    # this is the default
    --format genbank  
    # SeqIO format EMBL
    --format embl     

=back

=head1 Authors

Albert Vilella E<lt>avilella at ub dot eduE<gt>

=cut


use strict;

use Bio::DB::GenBank;
use Bio::SeqIO;
use Getopt::Long;

my $help = 0;            # WTH?
my $type = 'id';
my $outfmt = 'genbank';
my $featuretype = 'CDS';
my $onlygeneid = 0;

my $ok = GetOptions( 
		     'type:s' => \$type,
		     'outformat:s' => \$outfmt,
		     'h'           => \$help,
		     'help'        => \$help,
		     'onlygeneid'        => \$onlygeneid,
		     );

if((! $ok) || $help) {
    if(! $ok) {
	print STDERR "missing or unsupported option(s) on commandline\n";
    }
    system("perldoc $0");
    exit($ok ? 0 : 2);
}


my $gb = new Bio::DB::GenBank;

my $uid = shift @ARGV;

if ($type =~ /id/) {
    print STDERR "Getting file for $uid...\n";
    my $seqin = $gb->get_Seq_by_gi($uid);
    print "$uid;" , $seqin->display_id , ";" , $seqin->accession , ";";
    foreach my $f ( $seqin->get_SeqFeatures ) {
        if ($f->has_tag('db_xref')) {
            foreach my $tagval ($f->get_tag_values('db_xref')) {
                print "$tagval;";
            }
        }
    }
    print "\n";
    unless ($onlygeneid) {
        my $acc = $seqin->accession;
        my $output = $acc.".$outfmt";
    
        my $seqout = Bio::SeqIO->newFh ( -file   => ">$output",
                                         -format => $outfmt );
        print STDERR "File Output is: $output";
        print $seqout $seqin;
    }

} elsif ($type =~ /acc/) {
    print STDERR "Getting file for $uid...\n";
    my $seqin = $gb->get_Seq_by_acc($uid);
    my $acc = $seqin->accession;
    my $output = $acc.".$outfmt";
    
    my $seqout = Bio::SeqIO->newFh ( -file   => ">$output",
				     -format => $outfmt );
    print STDERR "File Output is: $output\n";
    print $seqout $seqin;

} elsif ($type =~ /batch/) {
    print STDERR "Getting file for $uid...\n";
    open (BATCH,"$uid") or die "$!";
    while (<BATCH>) {
        chomp $_;
        my $uid = $_;
        next unless ($uid =~ /\S+/);
        my $seqin = $gb->get_Seq_by_id($uid);
        print "$uid;" , $seqin->display_id , ";" , $seqin->accession , ";";
        foreach my $f ( $seqin->get_SeqFeatures ) {
            if ($f->has_tag('db_xref')) {
                foreach my $tagval ($f->get_tag_values('db_xref')) {
                    print "$tagval;";
                }
            }
        }
        print "\n";
        next if ($onlygeneid);
        my $acc = undef;
        eval {
            $acc = $seqin->accession;
            my $output = $acc.".$outfmt";
            my $seqout = Bio::SeqIO->newFh ( -file   => ">$output",
                                             -format => $outfmt );
            print STDERR "File Output is: $output\n";
            print STDERR "\n";
            print $seqout $seqin;
        };
        if ($@) {
            print STDERR "Error in $acc:\n";
            print STDERR "$@\n";
        }
    }
} else {
    die ("unknown --type of accession");
}



