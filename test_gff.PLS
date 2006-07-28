#!/usr/local/bin/perl

use Bio::DB::GFF;
use Getopt::Long;

my ($infile, $fastadir,$aggs,$output);
GetOptions(
	   'i|infile:s' => \$infile,
	   'f|fasta|fastadir:s'  => \$fastadir,
	   'a|agg|aggreg:s' => \$aggs,
	   'o|output:s'=> \$output,
          );

# my @aggs = split/\:/, $aggs || ('sequence','gene');

my $agg1 = Bio::DB::GFF::Aggregator->new
    (
     -source => 'assembly',
     -method => 'sequence',
    );
#                                          -sub_parts    => ['allele','variant']

my $agg2 = Bio::DB::GFF::Aggregator->new
    (
     -source       => 'pipeline',
     -sub_parts    => ['orthology_gene','orthology_transcript','orthology_exon','orthology_nonexon']
    );

# Open the sequence database
my $db      = Bio::DB::GFF->new
    (
     -adaptor => 'memory',
     -create  => 1,
     -aggregator => [$agg1,$agg2],
    );
#      '-aggregators' => \@aggs,

$db->initialize(1);
$db->load_gff("$infile");
# $db->load_fasta("$fastadir");

# # fetch a 1 megabase segment of sequence starting at landmark "ZK909"
# my $segment = $db->segment('ZK909', 1 => 1000000);

# fetch the chrX landmark
my $seg = $db->segment('4');
# my $seg = $db->segment('X', 1 => 100000);

# pull out all orthology_gene features
my @orthology_genes = $seg->features('orthology_gene');

1;;

# Get a region 500 bp upstream of the gene
for my $orth_gene (@orthology_genes) {
    my $upstream = $orth_gene->subseq(-1000,0);
    # get its DNA
    my $dna = $upstream->seq;
    # # and get all curated polymorphisms inside it
    # @polymorphisms = $upstream->contained_features('polymorphism:curated');
}

my @types = $db->types;

# count all feature types in the segment
my %type_counts = $segment->types(-enumerate=>1);

# get an iterator on all curated features of type 'exon' or 'intron'
my $iterator = $db->get_seq_stream(-type     => ['exon:curated','intron:curated']);

while (my $s = $iterator->next_seq) {
    print $s,"\n";
}

1;;

# # for each transcript, total the length of the introns
# my %totals;
# for my $t (@transcripts) {
#   my @introns = $t->Intron;
#   $totals{$t->name} += $_->length foreach @introns;
# }

# # Sort the exons of the first transcript by position
# my @exons = sort {$a->start <=> $b->start} $transcripts[0]->Exon;

# # Get a region 1000 bp upstream of first exon
# my $upstream = $exons[0]->subseq(-1000,0);

# # get its DNA
# my $dna = $upstream->seq;


