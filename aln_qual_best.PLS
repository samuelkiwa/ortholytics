#!/usr/local/bin/perl
use strict;
use Getopt::Long;

my ($inputfile,$debug);

GetOptions(
	   'i|input|inputfile:s' => \$inputfile,
           'debug' => \$debug,
          );

my $self = bless {};
open (INFILE, $inputfile) or die "$!\n";
while (<INFILE>) {
    $_ =~ /(\S+)\,\s+(\S+)\,\s+(\S+)\,\s+(\S+)\,\s+(\S+)\,\s+(\S+)\n/;
    next if $1 =~ /genes/;
    $self->{genes}{$1}{prog}{$2}{stats}{ovl} = $3;
    $self->{genes}{$1}{prog}{$2}{stats}{avg} = $4;
    $self->{genes}{$1}{prog}{$2}{stats}{pat} = $5;
    $self->{genes}{$1}{prog}{$2}{stats}{len} = $6;
}

my $bestprogavg = 'none';
my $bestprogovl = 'none';
my $bestprogpat = 'none';
my $bestovl = 0;
my $bestavg = 0;
my $bestpat = 0;
foreach my $gene (keys %{$self->{genes}}) {
    foreach my $prog (keys %{$self->{genes}{$gene}{prog}}) {
        if ($self->{genes}{$gene}{prog}{$prog}{stats}{avg} > $bestavg) {
            $bestprogavg = $prog;
            $bestavg = $self->{genes}{$gene}{prog}{$prog}{stats}{avg};
        }
        if ($self->{genes}{$gene}{prog}{$prog}{stats}{ovl} > $bestovl) {
            $bestprogovl = $prog;
            $bestovl = $self->{genes}{$gene}{prog}{$prog}{stats}{ovl};
        }
        if ($self->{genes}{$gene}{prog}{$prog}{stats}{pat} > $bestpat) {
            $bestprogpat = $prog;
            $bestpat = $self->{genes}{$gene}{prog}{$prog}{stats}{pat};
        }
        $self->{genes}{$gene}{best}{bestprogavg} = $bestprogavg;
        $self->{genes}{$gene}{best}{bestprogovl} = $bestprogovl;
        $self->{genes}{$gene}{best}{bestprogpat} = $bestprogpat;
        1;
    }
}

foreach my $gene (keys %{$self->{genes}}) {
    if ($self->{genes}{$gene}{best}{bestprogavg} eq $self->{genes}{$gene}{best}{bestprogovl}) {
        if ($self->{genes}{$gene}{best}{bestprogovl} eq $self->{genes}{$gene}{best}{bestprogpat}) {
            my $prog = $self->{genes}{$gene}{best}{bestprogavg};
            $self->{avg_ovl_pat}{$gene}{best}{$prog} = 1;
        }
    } elsif ($self->{genes}{$gene}{best}{bestprogavg} eq $self->{genes}{$gene}{best}{bestprogovl}) {
        my $prog = $self->{genes}{$gene}{best}{bestprogavg};
        $self->{avg_ovl_xxx}{$gene}{best}{$prog} = 1;
    } elsif ($self->{genes}{$gene}{best}{bestprogavg} eq $self->{genes}{$gene}{best}{bestprogpat}) {
        my $prog = $self->{genes}{$gene}{best}{bestprogavg};
        $self->{avg_xxx_pat}{$gene}{best}{$prog} = 1;
    } elsif ($self->{genes}{$gene}{best}{bestprogovl} eq $self->{genes}{$gene}{best}{bestprogpat}) {
        my $prog = $self->{genes}{$gene}{best}{bestprogovl};
        $self->{xxx_ovl_pat}{$gene}{best}{$prog} = 1;
    } elsif ($self->{genes}{$gene}{best}{bestprogavg} ne $self->{genes}{$gene}{best}{bestprogovl}) {
        if ($self->{genes}{$gene}{best}{bestprogovl} ne $self->{genes}{$gene}{best}{bestprogpat}) {
            $self->{xxx_xxx_xxx}{$gene}{best}{xxx} = 1;
        }
    }
}

print "avg_ovl_pat ", scalar(keys %{$self->{avg_ovl_pat}}), "\n";
print "avg_ovl_xxx ", scalar(keys %{$self->{avg_ovl_xxx}}), "\n";
print "avg_xxx_pat ", scalar(keys %{$self->{avg_xxx_pat}}), "\n";
print "xxx_ovl_pat ", scalar(keys %{$self->{xxx_ovl_pat}}), "\n";
print "xxx_xxx_xxx ", scalar(keys %{$self->{xxx_xxx_xxx}}), "\n";

1;
