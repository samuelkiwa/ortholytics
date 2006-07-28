# orthomcl_stats.PLS
#
# Cared for by Albert Vilella, Pablo Librado <>
#
# Copyright Albert Vilella, Pablo Librado
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

orthomcl_stats.PLS - DESCRIPTION 

=head1 SYNOPSIS

perl orthomcl_stats.PLS 
  -orth all_orthomcl.out
  -outfile myfile.csv
  -groupdir /my/dir
  [-thpercent 1]

=head1 DESCRIPTION

Describe the object here

=head1 AUTHOR - Albert Vilella, Pablo Librado

Email 

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut


# Let the code begin...

use strict;
use Getopt::Long;
use File::Path;

my ($orth,$outfile,$groupdir,$thpercent) = undef;

$thpercent = 1;

GetOptions(
	   'o|orth|ortho:s'  => \$orth,
	   'outfile:s'       => \$outfile,
	   'g|groupdir:s'       => \$groupdir,
	   'th|thpercent:s'       => \$thpercent,
          );

print "Reading orthomcl alignment sets...\n";

my $num_orth_rel = 0;
my $orthology;
my $regexp;
my $pattern;
my $numgenes;
my $numtaxa;
my $orth_string;

$pattern = '(ORTHOMCL(\S+)\((\d+)\ genes\,(\d+)\ taxa\)\:\s+(.+))';

$regexp = qr/^${pattern}/i;

my %ortho_rel;
my %names_taxa;

open (ORTH,"$orth") or die "cannot open orth file: $!";
print "Loading orthologies...\n";
while (<ORTH>) {
    if ($_ =~ /$regexp/g) {
        $orthology = sprintf("%05d", $2);
        $numgenes = $3;
        $numtaxa  = $4;
        $orth_string = $5;
        $ortho_rel{$orthology}{string} = $1;
        my @entries = split(/\ /, $orth_string);
        # Each ortho_rel will have the number of times a taxon is present
        foreach my $entry (@entries) {
            if ($entry =~ /\S+\((\S+)\)/) {
                $ortho_rel{$orthology}{species}{$1}++;
                $names_taxa{$1}++;
            }
#             else {
#                 warn "parsing error";
#             }
        }
        $orthology = "";
        $num_orth_rel++;
    }
}

close ORTH;

# Add a '0' if the taxon is missing in that ortho_rel
foreach my $orthology (keys %ortho_rel) {
    foreach my $taxa (sort keys %names_taxa) {
        unless (exists $ortho_rel{$orthology}{species}{$taxa}) {
            $ortho_rel{$orthology}{species}{$taxa} = 0;
        }
    }
}

# Orthology types strings - number of times each orthology type occurs
# Ex: 1-1-1-1- occurs X times
# Ex: 1-0-1-1- occurs Y times
print "Creating Orthology types strings...\n";
my $orth_type_string = undef;
my %types_orth;
foreach my $orthology (keys %ortho_rel) {
    foreach my $taxa (sort keys %names_taxa) {
        $orth_type_string .= "$ortho_rel{$orthology}{species}{$taxa}" . "-";
    }
    $types_orth{$orth_type_string}{occurrences}++;
    push @{$types_orth{$orth_type_string}{orthology}}, $orthology;
    $orth_type_string = '';
}
# Names taxa string
my $names_taxa_string = undef;
foreach my $taxa (sort keys %names_taxa) {
    $names_taxa_string .= "$taxa-";
}

open (OUTFILE,">$outfile") or die "cannot open orth file: $!";
print "Writing $outfile\n";


print OUTFILE "1 - Orthology types;\n";
print OUTFILE "$names_taxa_string;num_orth_rel;percent;\n";
my $percent = undef;
foreach my $type (sort { $types_orth{$b}{occurrences} <=> $types_orth{$a}{occurrences} } keys %types_orth) {
    $percent = sprintf ("%02.02f",($types_orth{$type}{occurrences}/$num_orth_rel)*100);
    print OUTFILE "$type;$types_orth{$type}{occurrences};$percent;\n";
    if ($thpercent <= $percent) {
        unless (-d "$groupdir/$type") {
            File::Path::mkpath("$groupdir/$type");
        }
        open (TYPE,">$groupdir/$type/$type.orth.txt") or die "cannot open file: $!";
            foreach my $orthology (@{$types_orth{$type}{orthology}}) {
                print TYPE "$ortho_rel{$orthology}{string}\n";
            }
            close TYPE;
    }
}
print OUTFILE "\n";

########################################

# Paralogy occurrences
my %types_para;
my $num_para_rel = 0;
my $has_inpara = 0;
foreach my $orthology (keys %ortho_rel) {
    foreach my $taxa (sort keys %names_taxa) {
        $orth_type_string .= "$ortho_rel{$orthology}{species}{$taxa}" . "-";
        # Keep cases were a taxon has more than 1 gene
        if (1 < $ortho_rel{$orthology}{species}{$taxa}) {
            $has_inpara = 1;
        }
    }
    if ($has_inpara) {
        my @para_occurrences = split(/\-/, $orth_type_string);
        my $num_nonORFans = 0;
        # Count the number of other taxa that is bigger than zero
        foreach my $occurrence (@para_occurrences) {
            if (0 < $occurrence) {
                $num_nonORFans++;
            }
        }
        # Keep as para if there are occurrences of other taxa:
        # Ex1: 2-1-0-0 this is a para
        # Ex1: 2-0-0-0 this is not a para
        if (1 < $num_nonORFans) {
            $types_para{$orth_type_string}++;
            $num_para_rel++;
        }
    }
    $orth_type_string = "";
    $has_inpara = 0;
}
print OUTFILE "2 - Only w/Paralogy Types;\n";
print OUTFILE "$names_taxa_string;num_para_rel;percent;\n";

$percent = undef;
foreach my $type (sort { $types_para{$b} <=> $types_para{$a} } keys %types_para) {
    $percent = sprintf ("%02.02f",($types_para{$type}/$num_para_rel)*100);
    print OUTFILE "$type;$types_para{$type};$percent;\n";
}

print OUTFILE "\n";

########################################

# ORFans types
print OUTFILE "3 - Only w/ORFans types;\n";

my %types_orfan;
my $num_orfan_rel = 0;
my $is_orfan = 0;
foreach my $orthology (keys %ortho_rel) {
    foreach my $taxa (sort keys %names_taxa) {
        $orth_type_string .= "$ortho_rel{$orthology}{species}{$taxa}" . "-";
        # Keep if there is any taxa missing
        if (0 == $ortho_rel{$orthology}{species}{$taxa}) {
            $is_orfan = 1;
        }
    }
    if ($is_orfan) {
        $types_orfan{$orth_type_string}++;
        $num_orfan_rel++;
    }
    $orth_type_string = "";
    $is_orfan = 0;
}
print OUTFILE "$names_taxa_string;num_orfan_rel;percent;\n";

$percent = undef;
foreach my $type (sort { $types_orfan{$b} <=> $types_orfan{$a} } keys %types_orfan) {
    $percent = sprintf ("%02.02f",($types_orfan{$type}/$num_orfan_rel)*100);
    print OUTFILE "$type;$types_orfan{$type};$percent;\n";
}

print OUTFILE "\n";

########################################

open (SUBCLUST,">$outfile.proposed_subclustering.txt") or die "cannot open proposed_subclustering file: $!";

print SUBCLUST "4 - Proposed subclusterings;\n";
print SUBCLUST "subclustering_string;num_sp;num_rel;\n";

my $subclustering_string;
# Order the ORFans by number of occurrences - higher occurrencies are
# first to propose for subclustering
foreach my $type (sort { $types_orfan{$b} <=> $types_orfan{$a} } keys %types_orfan) {
        my @occurrences = split(/\-/, $type);
        my @names_taxa  = split(/\-/, $names_taxa_string);
        my $pos;
        my $num_sp = 0;
        foreach my $sp (@occurrences) {
            unless (0 == $sp) {
                $subclustering_string .= "$names_taxa[$pos],";
                $num_sp++;
            }
            $pos++;
        }
        # Keep only when 3 or more species are present
        if (2 < $num_sp) {
            print SUBCLUST "$subclustering_string;$num_sp;$types_orfan{$type};\n";
        }
        $subclustering_string = "";
        1;
}

close SUBCLUST;
close OUTFILE;

print "--\n";
1;
