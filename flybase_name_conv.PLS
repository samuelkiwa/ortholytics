#!/usr/local/bin/perl

use Bio::Factory::FTLocationFactory;
use Getopt::Long;
my $fh;
my ($file,$gofiles,$all);

$all = 0;
GetOptions(
	   'i|infile:s' => \$file,
	   'g|gofiles:s'  => \$gofiles,
	   'all' => \$all,
          );

my @gofiles = split /\:/,$gofiles;

# open($fh, "grep '^>' $file") || die "$!";
open($fh, "$file") || die "$!";

my %entries;
my $id;

while (<$fh> ) {
    # >CG11023-PA
    if ( /^>(\S+)/ ) {
        $id = $1;
    }
    if ( /type=(\S+);/ ) {
        $entries{$id}{_type} = $1;
    }
    # loc=2L:join(7680..8116,8229..8589,8668..9273);
    if ( /loc=(\S+):(\S+);/ ) {
        my ($chrid,$locationstr) = ($1,$2);
        my $location = Bio::Factory::FTLocationFactory->from_string ($locationstr);
        my $runninglength = 0;
        my $i = 0;
        my @exons =  $location->each_Location;
        my $last = scalar @exons;
        for my $exon (@exons) {
            # I may be sloppy here, pls check that this is working the way you expect
            # defining A^TG is phase 1 and AT^G is phase 2 i
            my $phase = ( $runninglength += $exon->length) % 3;
            if ( $i != $last) {
#                 print "phase of intron $i is $phase\n";
            }
            $i++;
        }
    }
    # ID=CG11023-PA;
    if ( /ID=(\S+);/ ) {
        my $seqid = $1;
        $entries{$id}{_ID} = $1;
    }
    #  name=CG11023-PA;
    if ( /name=(\S+);/ ) {
    }
    #  db_xref=FlyBase:FBpp0088316,GB_protein:AAO41164.1,FlyBase:FBgn0031208,Gadfly:CG11023-PA;
    if ( /db_xref=(\S+);/ ) {
        my @xrefs = split /\,/,$1;
        foreach my $xref (@xrefs) {
            my ($key,$value) = split /\:/,$xref;
            $entries{$id}{_xref}{$key} = $value;
        }
    }
    #  species=dmel;
    if ( /species=(\S+)/ ) {
        $entries{$id}{_species} = $1;
    }
    #  len=1404
    if ( /len=(\S+)/ ) {
        $entries{$id}{_len} = $1;
    }
}

my %flybase;

foreach my $entry (keys %entries) {
    my $fbgn = $entries{"$entry"}{'_xref'}{'FlyBase'};
    if ($all) {
        $flybase{$fbgn} .= "$entry:";
    } else {
        if ($entry =~ /-PA/) {
            $entry =~ s/-PA//g;
            $flybase{$fbgn} .= "$entry" 
        }
    }
}

foreach my $gofile (@gofiles) {
    open (FILE, "$gofile") or die "$!";
    my %seen =();
    while (<FILE>) {
#         while ( /(\w[\w'-]*)/g ) {
        foreach my $word (split) {
            $seen{$word}++;
            if ($word =~ /GO\:(\S+)/) {
                my $filename = "$gofile.GO"."$1";
                open (NEW, ">$filename") or die "$!";
                print NEW "#$word\n#" if ($word =~ /GO\:/);
            } elsif (!defined $flybase{$word}) {
                print NEW "#$word\n";
            } else {
                $word =~ s/(\S+)/$flybase{$word}/g if (defined $flybase{$word});
                print NEW "$word\n";
            }
        }
        close(NEW);
    }
    close(FILE);
    foreach my $word (sort { $seen{$b} <=> $seen{$a} } keys %seen) {
        printf "%5d %s\n", $seen{$word}, $word;
    }
}

1;;
