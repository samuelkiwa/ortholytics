#!/usr/local/bin/perl

use strict;
use Getopt::Long;
use GO::AppHandle;

my ($alndir,$file,$dbname,$mysqlhost,$all);

GetOptions(
	   'alndir:s' => \$alndir,
	   'i|infile:s' => \$file,
	   'd|database:s' => \$dbname,
	   'h|host:s' => \$mysqlhost,
	   'all' => \$all,
          );

my @alnlist;
opendir ALNDIR, $alndir or die "couldnt find $alndir:$!";
foreach my $name (readdir ALNDIR) {
    next if $name =~ /^\.\.?$/; # skip over dot files
    next unless -d "$alndir/$name";
    next unless $name =~ /CG/;
    push @alnlist, $name;
  }
closedir ALNDIR;

my $fh;
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
#     if ( /loc=(\S+):(\S+);/ ) {
#         my ($chrid,$locationstr) = ($1,$2);
#         my $location = Bio::Factory::FTLocationFactory->from_string ($locationstr);
#         my $runninglength = 0;
#         my $i = 0;
#         my @exons =  $location->each_Location;
#         my $last = scalar @exons;
#         for my $exon (@exons) {
#             # I may be sloppy here, pls check that this is working the way you expect
#             # defining A^TG is phase 1 and AT^G is phase 2 i
#             my $phase = ( $runninglength += $exon->length) % 3;
#             if ( $i != $last) {
# #                 print "phase of intron $i is $phase\n";
#             }
#             $i++;
#         }
#     }
    # ID=CG11023-PA;
    if ( /ID=(\S+);/ ) {
        my $seqid = $1;
        $entries{$id}{_ID} = $1;
    }
    #  name=CG11023-PA;
    if ( /name=(\S+);/ ) {
        $entries{$id}{_name} = $1;
    }
    #  db_xref=FlyBase:FBpp0088316,GB_protein:AAO41164.1,FlyBase:FBgn0031208,Gadfly:CG11023-PA;
    if ( /db_xref=(\S+);/ ) {
        my @xrefs = split /\,/,$1;
        foreach my $xref (@xrefs) {
            my ($key,$value) = split /\:/,$xref;
            $entries{$id}{_xref}{$key} = lc($value);
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
    $fbgn = lc($fbgn);
    if ($all) {
        $flybase{$entry} = "$fbgn";
    } else {
        if ($entry =~ /-PA/) {
            $entry =~ s/-PA//g;
            if ($fbgn =~ /fbgn/i) {
                $fbgn =~ s/fbgn/FBgn/;
                $flybase{$entry} = "$fbgn";
            }
        }
    }
}

my %aln_present_flybase;
foreach my $aln (@alnlist) {
    unless ($all) {
        $aln =~ s/-P\S//g;
    }
    if (defined $flybase{$aln}) {
        $aln_present_flybase{$aln} = $flybase{$aln};
    }
}

my @glist = values %aln_present_flybase;

# connect to a database on a specific host
my $apph = GO::AppHandle->connect
    (
     -dbname=>$dbname,
     -dbhost=>$mysqlhost
    );

# EXAMPLE 5
# fetching a subgraph of GO,
# constraining by gene products

# get all terms that were used to annotate these two SGD genes
# my $graph = $apph->get_graph_by_search("galectin", -1);

#      {products=>["Eip63F-1", "Krt1-13"]}
#      {products=>\@glist}

my $products = $apph->get_product
    (
     {
      accs=>\@glist
     }
    );

my $terms = 
    $apph->get_terms
    (
     {
      products=>[$products]
     }
    );

# my $term_l = 
#     $apph->get_terms
#     (
#      {
#       products=>\@glist,
#       search_fields=>"name,synonym,dbxrefs,comments",
#      }
#     );

#       search_fields=>"name,synonym,dbxrefs,comments",


# build a graph all the way to the leaf nodes
# from the above terms
my $graph = $apph->get_graph_by_terms
    (
     $terms,
     -1
    );

print @$graph;

exit;

1;;

# create an iterator on the graph
my $it = $graph->create_iterator;

# iterate through every node in graph
while (my $ni = $it->next_node_instance) {
    my $depth = $ni->depth;
    my $term = $ni->term;
    printf 
    "%s Term = %s (%s)  // ASSOCS=%s\n",
    "----" x $depth,
    $term->name,
    $term->public_acc,
    join("; ",
    map {$_->gene_product->acc} @{$term->association_list});
}

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

# while (num_nodes(tree) > 1) {
#     level = biggest(graph);
#     # 1) trim the lowest that >min
#     for (nodes in tree(level)) {
#         if ((no has_child) && (>min) ) {
#             trim_save;
#             trim_from_parents;
#         }
#     }
#     # 2) merge to parent the lowest that <=min and only one parent
#     for (nodes in tree(level)) {
#         if ((no has_child) && ('<='min) && (1 == num_parents(node)) ) {
#             trim;
#         }
#     }
#     # FIXME: decide wisely in diagram1 and diagram2
#     # 3a) merge to parent the lowest that <=min and more than one parent
#     for (nodes in tree(level)) {
#         if ((no has_child) && ('<='min) && (1 < num_parents(node)) {
#             if (merge to smallest parent that after merging reaches >min) {
#                 merge to smallest parent that after merging reaches >min;
#             }
#         }
#     }
#     # 3b)
#     for (nodes in tree(level)) {
#         if ((no has_child) && ('<='min) && (1 < num_parents(node)) {
#             merge to smallest parent;
#         }
#     }
# }
# trim_save($node_on_top)



# EXAMPLE 1
# fetching a GO term from the datasource

# my $term = $apph->get_term({acc=>"GO:0003677"});
# printf 
#     "GO term; name=%s GO ID=%s\n",
#     $term->name(), $term->public_acc();

# EXAMPLE 2
# fetching a list of associations to the ER
# (and all the GO terms that are subtypes of the ER, or 
#  located within the ER)
# for which there is reasonably good evidence
# (traceable author / direct assay)

# my $assocs = $apph->get_associations
#     ({name=>"endoplasmic reticulum"},
#      {
#       evcodes=>["TAS", "IDA"]});
# foreach my $assoc (@$assocs) {
#     printf
#         "Gene: %s evidence for association: %s %s",
#             $assoc->gene_product->symbol,
#                 $assoc->evidence->code(),
#                     $assoc->evidence->xref->xref_key();
# }
# EXAMPLE 3
# fetching a subgraph of GO

# my $graph = $apph->get_graph
#     (
#      -acc=>3677,
#      -depth=>3
#     );
# foreach my $term (@ {$graph->get_all_nodes}) {
#     printf 
#         "GO term; name=%s GO ID=%s\n",
#             $term->name(), $term->public_acc();
# }

# EXAMPLE 4
# fetching a subgraph of GO,
# and using a graph iterator to 
# display the graph

# my $graph = $apph->get_graph_by_search
#     (
#      "DNA helicase*"
#     );
# my $it = $graph->create_iterator;

# while (my $ni = $it->next_node_instance) {
#     my $depth = $ni->depth;
#     $term = $ni->term;
#     printf 
#         "%s Term = %s (%s)  // n_assocs=%s // depth=%d\n",
#             "----" x $depth,
#                 $term->name,
#                     $term->public_acc,
#                         $term->n_associations || 0,
#                             $depth;
# }


