#!/usr/bin/perl

# ancestors_go.PLS
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

ancestors_go.PLS - DESCRIPTION 

=head1 SYNOPSIS

perl ancestors_go.PLS  \
    -i /home/avb/wallace/eukarya/drosophila/funccats/flybase/GO.txt \

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
    $godir,
    $plain,
   );

$plain = 0;
$godir = "/home/avb/wallace/go";

GetOptions(
	   'i|input|inputfile:s' => \$inputfile,
	   'g|godir:s' => \$godir,
	   'l|list|genelist:s' => \$genelist,
           'plain' => \$plain,
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
print "# Manually loading annotation file...\n";
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

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst);
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $start_time = sprintf ("%04d%02d%02d_%02d%02d%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec);
print "# $start_time\n";

print "# Loading ontology (might take a couple of minutes)...\n";
$parser = new GO::Parser({handler=>'obj'});
if ($plain) {$parser->parse("$godir/gene_ontology.obo");}
else {$parser->parse("$godir/go_daily-termdb.obo-xml");}
$graph = $parser->handler->graph;

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $end_time = sprintf ("%04d%02d%02d_%02d%02d%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec);
print "# $end_time\n";

# print "# Loading association (might take a couple of minutes)...\n";
# $parser->parse("$inputfile"); # gene assocs
# ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
# my $end_time = sprintf ("%04d%02d%02d_%02d%02d%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec);
# print "# $end_time\n";

my %groups;
my %leaf_groups;
foreach my $entry (keys %entries) {
    next unless ($entry =~ /GO/);
    1;;
    $term = $graph->get_term($entry);
    if (!defined $term) {
        next; # do nothing right now
#         # FIXME: Give me the base term for this alt_id
        print "# error $entry\n";
#         my $base_id = $graph->base_id($entry); # get the base_id here
#         print "#   base $base_id\n";
#         $term = $graph->get_term($base_id);
    }
    # fill groups hash
    my $term_acc = $term->public_acc;
    my $term_name = $term->name;
    $leaf_groups{$term_acc}{_totalcount}++;
    $leaf_groups{$term_acc}{_name} = $term_name;
    $leaf_groups{$term_acc}{_level} = 0;
    $ancestor_terms =
        $graph->get_recursive_parent_terms($term->public_acc);
    my $rel_level = 1;
    foreach my $anc_term (@$ancestor_terms) {
        my $anc_term_acc = $anc_term->public_acc;
        $leaf_groups{$term_acc}{_road} .= "$anc_term_acc:";
        my $anc_term_name = $anc_term->name;
#         my $n_deep_assoc = $graph->n_deep_associations($anc_term_acc);
#         print "# deep_assoc $term_acc , $n_deep_assoc\n";
        $groups{$anc_term_acc}{_totalcount}++;
        $groups{$anc_term_acc}{_name} = $anc_term_name;
#         if ($groups{$anc_term_acc}{_minlevel} > $rel_level && $rel_level > 0) {
#             $groups{$anc_term_acc}{_minlevel} = $rel_level;
#             $groups{$anc_term_acc}{_minlevelid} = $term_acc;
#         }
        if ($groups{$anc_term_acc}{_maxlevel} < $rel_level) {
            $groups{$anc_term_acc}{_maxlevel} = $rel_level;
            $groups{$anc_term_acc}{_maxlevelid} = $term_acc;
        }
        $rel_level++;
    }
}

# # This is only for pretty-print debugging
# print "#totalcount ; maxlevel ; maxlevelid ; group ; name\n";
# my @entries;
# foreach my $group (keys %groups) {
#     1;;
#     my $totalcount = sprintf ("%05d", $groups{$group}{_totalcount});
#     my $name = $groups{$group}{_name};
#     my $maxlevel = sprintf ("%03d", $groups{$group}{_maxlevel});
#     my $maxlevelid = sprintf ("%03d", $groups{$group}{_maxlevelid});
#     push @entries, "$totalcount ; $maxlevel ; $maxlevelid ; $group ; $name\n";
# }

# print sort {$b <=> $a} @entries;

# rules (most important to less important)

# a) >min
# b) <max

# i) closest to median to decide which parent (in case of multiple
# options)
# ii) chosen nodes are the lowest possible


# while (E nodes in tree)
# 1) trim the lowest that >min
# 2) merge to parent the lowest that <=min and only one parent
# 3) merge to parent the lowest that <=min
     # FIXME: recurse over parents to decide - very long
#    a) merge to smallest parent that after merging reaches >min
#       if tie (random)
#    b) merge to smallest parent
# -- loop

while (num_nodes(tree) > 1) {
    level = biggest(graph);
    # 1) trim the lowest that >min
    for (nodes in tree(level)) {
        if ((no has_child) && (>min) ) {
            trim_save;
            trim_from_parents;
        }
    }
    # 2) merge to parent the lowest that <=min and only one parent
    for (nodes in tree(level)) {
        if ((no has_child) && ('<='min) && (1 == num_parents(node)) ) {
            trim;
        }
    }
    # FIXME: decide wisely in diagram1 and diagram2
    # 3a) merge to parent the lowest that <=min and more than one parent
    for (nodes in tree(level)) {
        if ((no has_child) && ('<='min) && (1 < num_parents(node)) {
            if (merge to smallest parent that after merging reaches >min) {
                merge to smallest parent that after merging reaches >min;
            }
        }
    }
    # 3b)
    for (nodes in tree(level)) {
        if ((no has_child) && ('<='min) && (1 < num_parents(node)) {
            merge to smallest parent;
        }
    }
}
trim_save($node_on_top)




foreach my $group (keys %leaf_groups) {
    print "$group = $leaf_groups{$group}{_road}\n";
}

1;;


