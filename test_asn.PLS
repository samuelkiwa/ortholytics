#!/usr/local/bin/perl
use Getopt::Long;
use Data::Dumper;

my ($inputfile,$yesnoopt);

GetOptions(
	   'i|input|inputfile:s' => \$inputfile,
           'yesnoopt' => \$yesnoopt,
          );

my $parser = Bio::ASN1::Sequence->new('file' => $inputfile);
1;;
while (my $result = $parser->next_seq) {
    1;;
    Dumpvalue->new->dumpValue($result);
}


# a new way to get the $result data hash for a particular sequence id:
use Bio::ASN1::Sequence::Indexer;
my $inx = Bio::ASN1::Sequence::Indexer->new(-filename => 'seq.idx');
my $seq = $inx->fetch_hash('AF093062');

# for creation of .idx index files please refer to
# Bio::ASN1::Sequence::Indexer perldoc
