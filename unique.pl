#!/usr/local/bin/perl
foreach ( <> ) { $unique{$_} = 1 ; }
print keys(%unique); # values(%unique) is the other half
