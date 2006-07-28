use strict;
use Bio::SeqIO;

my $file = shift;

my $seqio = Bio::SeqIO->new(-format => 'entrezgene',
                            -file => $file);
my $gene = $seqio->next_seq;
