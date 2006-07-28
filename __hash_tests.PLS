#!/usr/local/bin/perl

# foreach ( <> ) { $unique{$_} = 1 ; }
# print keys(%unique); # values(%unique) is the other half
########################################

my @A = (1,2,2,3,3,6,5,4,3,2,23,23,4,6,6,5,2,234,123,4124,21,43,124,1234);
my @B = (1,21,4,43,1,121,331,34,4,114,3321,1,43,34,44,3,3,3,3,2,2,3,233,32,32);
my %seen; # lookup table
my @aonly;# answer
my @bonly;# answer

my %bseen;
foreach $item (@B) {
    $bseen{$item} = 1;
}

# build lookup table
@seen{@B} = ();

foreach $item (@A) {
    push(@aonly, $item) unless exists $seen{$item};
    push(@bonly, $item) unless exists $bseen{$item};
}

1;
