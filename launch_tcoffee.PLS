#!/usr/local/bin/perl
# launch_muscle.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

launch_muscle.PLS -  

=head1 SYNOPSIS

perl launch_muscle.PLS -dir /my/dir

=head1 DESCRIPTION

launch_muscle will load the sequences for all the one-seq files in the
specified format in the speficied

=head1 AUTHOR - Albert Vilella

Email 

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut


# Let the code begin...

use strict;
use Bio::Tools::Run::Alignment::Muscle;
use Bio::SeqIO;

# Build a muscle alignment factory
@params = ('ktuple' => 2, 'matrix' => 'BLOSUM');
$factory = new Bio::Tools::Run::Alignment::Muscle (@params);

while (defined($file = readdir($dir))) {
    if ($file =~ /(matrix\S+(NC\S+)\.txt)$/) {
        $matrix_names{$2} = $1;
    } elsif ($file =~ /((NC\S+)\.habitat\.txt)$/) {
        $habitat_names{$2} = $1;
    }
}

eval {
    my $input = Bio::SeqIO->new
        (-file=>Bio::Root::IO->catfile("t","data","test.waba"));
    ok my $seq = $input->next_seq();
};
$@ ? ok 1 : ok 0;


my $muscle_present = $factory->executable();
unless ($muscle_present) {
    warn "muscle program not found. Skipping tests $Test::ntest to $NUMTESTS.\n";
    exit(0);
}

# Pass the factory an array of Bio::Seq objects
$seq_array_ref = \@seq_array;
$aln = $factory->align($seq_array_ref);



1;
