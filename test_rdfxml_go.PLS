# use GO::Parser;
use GO::AppHandle;
use GO::IO::RDFXML;

my $apph = GO::AppHandle->connect(-d=>$go, -dbhost=>$dbhost);
my $term = $apph->get_term({acc=>00003677});

#### ">-" is STDOUT
my $out = new FileHandle(">-");  

my $xml_out = GO::IO::XML->new(-output=>$out);
$xml_out->start_document();
$xml_out->draw_term($term);
$xml_out->end_document();
