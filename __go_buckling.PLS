#!/usr/local/bin/perl

use IO::File;
use POSIX qw(tmpnam);

$MappedAnnotFile = MapToParents($AnnotationFile);
my (%TermFreq, %TermName, %TermCateg) = ParseAnnotation($MappedAnnotFile);

foreach my $gene (@Genelist) {
    my ($annot, $Cat) = GetAnnotation($gene);
    my %chosenSlimTerms;
    $chosenSlimTerms{$annot} = "";
    %chosenSlimTerms = DecideSlim($annot, $Cat,$acceptableidsAtTerm, \%chosenSlimTerms);
    # The %chosenSlimTerms has the fine annotation as its key and
    # the tab delimited list of chosen slim terms as the value
    #         ... Do something with the $gene and %chosenslimterms ...
    #         ... May be write to a file called "custom annotations" and analyze it...
    #         ... You can visualize the "cut" in the DAG by painting the fine annotation
    #             in one color and its chosen slim parent in another color.
}

sub MapToParents{
    my $FullAnnotFile = $_[0];
    my ($FH2,$PAnnotFile) = tmpnam();

    my $dbh = DBI->connect($data_source,$username,$password) or
        die ("Can't connect to $data_source:$@");

    open (File, $FullAnnotFile) or die "cant open $FullAnnotFile";
    open ($FH2, ">$PAnnotFile") or die "cant write output";
    my %Unique;
    while (<File>) {
        chomp $_;
        my ($Mips, $GeneModel, $GO_id, $GO_annot, $undef, $GO_categ) = split(/\t/,$_);
        next if ($GO_id=~m/null/i);
        # Find the Parents of a term -- (includes the term too):
        my ($ParentIds,$ParentAnnot)=GetGoParents($GO_id,$GO_categ,$dbh);
        for my $i (0..$#{$ParentIds}) {
            if (!defined $Unique{lc($Mips.$$ParentIds[$i])}) {
                print $FH2 "$Mips\t$GeneModel\t$$ParentIds[$i]\t"
                    ."$$ParentAnnot[$i]\t$undef\t$GO_categ\n";
                $Unique{lc($Mips.$$ParentIds[$i])}=1;
            }
        }
    }
    close File;
    close $FH2;
    undef %Unique;
    $dbh->disconnect;
    return $PAnnotFile;
}

sub ParseAnnotation {
    # Parse your fine level annotations file (tab delmited text similar to the 
    # go_associations file) to find the number of genes at each term

    my $AnnotFile = $_[0];
    open (File, $AnnotFile) or die "cant open $AnnotFile";
    my (%TermFreq, %TermName, %TermCateg);
    while (<File>) {
        chomp $_;
        my ($gene, $GeneModel, $GO_id, $GO_annot, $undef, $GO_categ) = split(/\t/,$_);
        $gene =~s/^\s+//;
        if (defined $TermFreq{$GO_id}) {
            $TermFreq{$GO_id}++;
        } else {
            $TermFreq{$GO_id}=1;
        }
        $TermName{$GO_id} = $GO_annot;
        $TermCateg{$GO_id} = $GO_categ;
    }
    close File;
    return (\%TermFreq, \%TermName, \%TermCateg);
}

sub GetAnnotation {}
#  -- will just fetch the finest level annotation and its category for $gene

sub DecideSlim {
    # recursively keeps calling itself till it finds a suitable parent
    # that has an acceptable number of genes assigned to it.

    my ($fineannot, $Cat, $cutoff, $chosenSlimTermsRef) = @_;
    my ($parentsRef, $parents_labelRef) = GetParent($fineannot,$Cat,$dbh);
    for my $i (0..$#{$parentsRef}) {
        my $p = $$parentsRef[$i];
        # I am using the %TermFreq created by sub ParseAnnotation
        if ($$TermFreq{$p} >= $cutoff) {
            if ($chosenSlimTerms{$fineannot}) {
                $$chosenSlimTermsRef{$fineannot}= $$TermFreq{$p};
            } else {
                $$chosenSlimTermsRef{$fineannot} .= "\t".$$TermFreq{$p};
            }
        } else {
            %chosenSlimTerms = DecideSlim ($p, $Cat, $cutoff, $chosenSlimTermsRef);
        }
    }
    return %chosenSlimTerms;
}

sub GetGOParents {
    my ($GO_id,$GO_categ,$dbh)=@_;
    my (@Parent_annot,@Parent_id,$GO_categ_full);

    my $query = "select p.acc, p.name from graph_path INNER JOIN term AS t ON "
        ."(t.id = graph_path.term2_id) INNER JOIN term AS p ON "
            ."(p.id = graph_path.term1_id) where t.acc = ? and "
                ."graph_path.distance =1 and t.term_type = ? ";

    if ($GO_categ =~m/process/i) {
        $GO_categ_full = "biological_process";
    }
    if ($GO_categ =~m/component/i) {
        $GO_categ_full = "cellular_component";
    }
    if ($GO_categ =~m/function/i) {
        $GO_categ_full = "molecular_function";
    }

    my $sth = $dbh->prepare("$query");
    $sth->execute($GO_id, $GO_categ_full);
    while (my($GO_ids, $GO_annots) = $sth->fetchrow_array) { 
        push (@Parent_id, $GO_ids);
        push (@Parent_annot, $GO_annots);
    }
    return (\@Parent_id,\@Parent_annot);
}

