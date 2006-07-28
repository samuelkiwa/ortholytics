#!/usr/local/bin/perl
push @a, 1..100;

@b = ();
for (@a) {
    my $i = rand @b+1;
    push @b, $b[$i];
    $b[$i] = $_;
}

print "\n";
1;
