# FETCHING GRAPH FROM FILES
use GO::Parser;
my $parser = new GO::Parser({handler=>'obj'});
$parser->parse("/home/avb/wallace/eukarya/drosophila/funccats/flybase/GO.txt"); # gene assocs
$parser->parse("/home/avb/wallace/go/gene_ontology.obo");     # ontology
# get L<GO::Model::Graph> object
my $graph = $parser->handler->graph;
my $terms = $graph->term_query("/kinase/");  # matching terms
foreach my $term (@$terms) {
  # find gene products associated to this term
  my $assocs = $graph->deep_association_list($term->acc);
  printf "Term: %s %s\n", $term->acc, $term->name;
  print "  Associations (direct and via transitive closure_\n";
  foreach my $assoc (@$assocs) {
    next if $assoc->is_not;
    printf "  Assoc evidence: %s to: %s %s\n",
      join(';', map {$_->code} @{$assoc->evidence_list}),
      $assoc->gene_product->xref->as_str,
      $assoc->gene_product->symbol;
  }
}
