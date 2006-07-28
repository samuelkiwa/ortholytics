#!/usr/bin/perl

# create_go_genelists.PLS
#
# Cared for by Albert Vilella <>
#
# Based on specification of the gene_association format is defined at:
#   http://www.geneontology.org/GO.annotation.html#file
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

create_go_genelists.PLS - DESCRIPTION 

=head1 SYNOPSIS

perl create_go_genelists.PLS  \
    -i /home/avb/wallace/eukarya/drosophila/funccats/flybase/GO.txt \
    -outdir /home/avb/wallace/eukarya/drosophila/funccats/flybase/genegroups \

=head1 DESCRIPTION

Describe the object here

=head1 AUTHOR - Albert Vilella

Email 

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut

# Let the code begin...

use strict;
use Getopt::Long;
use GO::Parser;
use File::Path;

my ($inputfile,
    $genelist,
    $outdir,
    $godir,
   );

$godir = "/home/avb/wallace/go";

GetOptions(
	   'i|input|inputfile:s' => \$inputfile,
	   'o|outdir:s' => \$outdir,
	   'g|godir:s' => \$godir,
	   'l|list|genelist:s' => \$genelist,
          );

my $linenum = 0;
my $line = "";
my @errors = ();
my @column = ();

#
# Column positions 
#
use constant DB => 0;
use constant DB_OBJECT_ID => 1;
use constant DB_OBJECT_SYMBOL => 2;
use constant QUALIFIER => 3;
use constant GOID => 4;
use constant REFERENCE => 5;
use constant EVIDENCE => 6;
use constant WITH => 7;
use constant ASPECT => 8;
use constant DB_OBJECT_NAME => 9;
use constant DB_OBJECT_SYNONYM => 10;
use constant DB_OBJECT_TYPE => 11;
use constant TAXON => 12;
use constant DATE => 13;
use constant ASSIGNED_BY => 14;

# Number of TAB delimited columns in file
use constant COLNUM => 15;

$column[DB] = ["DB", 1, 1];
$column[DB_OBJECT_ID] = ["DB_Object_ID", 1, 0];
$column[DB_OBJECT_SYMBOL] = ["DB_Object_Symbol", 1, 0];
$column[QUALIFIER] = ["Qualifier", 0, 1];
$column[GOID] = ["GOID", 1, 1];
$column[REFERENCE] = ["DB:Reference", 1, 1];
$column[EVIDENCE] = ["Evidence", 1, 1];
$column[WITH] = ["With", 0, 0];
$column[ASPECT] = ["Aspect", 1, 0];
$column[DB_OBJECT_NAME] = ["DB_Object_Name", 0, 0];
$column[DB_OBJECT_SYNONYM] = ["DB_Object_Synonym", 0, 0];
$column[DB_OBJECT_TYPE] = ["DB_Object_Type", 1, 1];
$column[TAXON] = ["Taxon", 1, 1];
$column[DATE] = ["Date", 1, 1];
$column[ASSIGNED_BY] = ["Assigned_by", 1, 0];

my (@cols,%entries);

open(INPUTFILE,"$inputfile") or die "couldnt open $inputfile: $!";
print "Loading annotation file...";
while ( defined($line = <INPUTFILE>) ) {
    $linenum++;
    unless ( $line =~ /.*\n/ ) {
        die "error in file: $!";
    }
	
    chomp $line;

# skip comment lines
    next if ($line =~ m/^\!/);

# blank line?
    if ( $line eq "" ) {
        die "error in file: $!";
	next;
    }

# split TAB delimited columns
    @cols = split(/\t/, $line);

    if ( @cols ne COLNUM) {
        die "invalid number of columns in file: $!";
    }

    if ($cols[GOID]) {
        $entries{$cols[GOID]}{_DB} = $cols[DB];
        $entries{$cols[GOID]}{_DB_OBJECT_ID} = $cols[DB_OBJECT_ID];
        $entries{$cols[GOID]}{_DB_OBJECT_SYMBOL} = $cols[DB_OBJECT_SYMBOL];
        $entries{$cols[GOID]}{_QUALIFIER} = $cols[QUALIFIER];
        $entries{$cols[GOID]}{_GOID} = $cols[GOID];
        $entries{$cols[GOID]}{_REFERENCE} = $cols[REFERENCE];
        $entries{$cols[GOID]}{_EVIDENCE} = $cols[EVIDENCE];
        $entries{$cols[GOID]}{_WITH} = $cols[WITH];
        $entries{$cols[GOID]}{_ASPECT} = $cols[ASPECT];
        $entries{$cols[GOID]}{_DB_OBJECT_NAME} = $cols[DB_OBJECT_NAME];
        $entries{$cols[GOID]}{_DB_OBJECT_SYNONYM} = $cols[DB_OBJECT_SYNONYM];
        $entries{$cols[GOID]}{_DB_OBJECT_TYPE} = $cols[DB_OBJECT_TYPE];
        $entries{$cols[GOID]}{_TAXON} = $cols[TAXON];
        $entries{$cols[GOID]}{_DATE} = $cols[DATE];
        $entries{$cols[GOID]}{_ASSIGNED_BY} = $cols[ASSIGNED_BY];
    }
}

close INPUTFILE;

my ($parser,$graph,$term,$ancestor_terms);
print "Loading ontology (might take a couple of minutes)...";
$parser = new GO::Parser({handler=>'obj'});
$parser->parse("$godir/gene_ontology.obo");
$graph = $parser->handler->graph;

unless (-d "$outdir") {
    File::Path::mkpath("$outdir");
}

my %groups;
my ($term_acc,$term_name,$term_string);
foreach my $entry (keys %entries) {
    next unless ($entry =~ /GO/);
    $term = $graph->get_term("$entry");
    $graph = $parser->handler->graph;
    printf "# Term: %s;%s\n", $term->acc, $term->name;
    $term_acc = $term->acc;$term_name = $term->name;
    $term_string = "# Term: $term_acc;$term_name\n";
    $ancestor_terms =
        $graph->get_recursive_parent_terms($term->acc);
    foreach my $anc_term (@$ancestor_terms) {
#         printf "  Ancestor term: %s %s\n", $anc_term->acc, $anc_term->name;
        my $anc_term_acc = $anc_term->acc;
        my $num = sprintf("%05d",1);
        if($entries{$anc_term_acc}{_DB_OBJECT_SYMBOL}) {
            my $symbol = $entries{$anc_term_acc}{_DB_OBJECT_SYMBOL};
            unless ($symbol =~ /CG\d+/) {
                my @syns = split(/\|/, $entries{$anc_term_acc}{_DB_OBJECT_SYNONYM});
                foreach my $syn (@syns) {
                    if ($syn =~ /CG\d+/) {
                        $symbol = $syn;
                    }
                }
            }
            $groups{$anc_term_acc}{$num} = $symbol;
            $num++; $num = sprintf("%05d",$num);
        } else {
            #FIXME: synonyms
        }
    }
    my $no_bashlashed_entry = $entry;
    $no_bashlashed_entry =~ s/\://;
    open(OUTFILE,">$outdir/$no_bashlashed_entry") or die "couldnt create file: $!\n";
    print OUTFILE "$term_string";
     foreach my $group (keys %groups) {
        foreach my $entry (keys %{ $groups{$group} }) {
            print OUTFILE "$groups{$group}{$entry}\n";
        }
    }
    close OUTFILE;
}

1;;


