#!/usr/bin/perl -w
use strict;
use Bio::Tools::Phylo::PAML;

my $outcodeml = shift(@ARGV);

my $paml_parser = new Bio::Tools::Phylo::PAML(-file => "./$outcodeml", 
					      -dir => "./");
my $result = $paml_parser->next_result();
my $MLmatrix = $result->get_MLmatrix(); # get MaxLikelihood Matrix
my @otus = $result->get_seqs;
if( $#{$MLmatrix} < 0 ) {
    for my $tree ($result->next_tree ) {
	for my $node ( $tree->get_nodes ) {
	    my $id;
	    if( $node->is_Leaf() ) {
		$id = $node->id;
	    } else {
		$id = "(".join(",", map { $_->id } grep { $_->is_Leaf } 
			   $node->get_all_Descendents) .")";
	    }
	    if( ! $node->ancestor || ! $node->has_tag('t') ) {
		# skip when no values have been associated with this node
		# (like the root node)
		next; 
	    }
            # I know this looks complicated
	    # but we use the get_tag_values method to pull out the annotations
	    # for each branch
	    # The ()[0] around the call is because get_tag_values returns a list
	    # if we want to just get the 1st item in the list we have 
	    # to tell Perl we are treating it like an array.
	    # in the future get_tag_values needs to be smart and just
	    # return the 1st item in the array if called in scalar
	    # context
	    
	    printf "%s\tt=%.3f\tS=%.1f\tN=%.1f\tdN/dS=%.4f\tdN=%.4f\tdS=%.4f\tS*dS=%.1f\tN*dN=%.1f\n",
	    $id,
	    map { ($node->get_tag_values($_))[0] }
	    qw(t S N dN/dS dN dS), 'S*dS', 'N*dN';
	}
    }
} else {
    my $i =0;
    my @seqs = $result->get_seqs;
    for my $row ( @$MLmatrix ) {
	print $seqs[$i++]->display_id, join("\t",@$row), "\n";
    }
}
