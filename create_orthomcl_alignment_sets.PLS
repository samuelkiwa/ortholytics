# create_orthomcl_alignment_sets.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

create_orthomcl_alignment_sets.PLS - open a list of
1-1(-1-1-1-1-1-1-1-1-1-1) orthologies and the multifasta nucleotide
files with the ">locus_id_" description, and create a directory with
each set of fastas for each entry in the orthologies list.

=head1 SYNOPSIS

create_orthomcl_alignment_sets.PLS
 -i1 file1.fna
 -i2 file2.fna
 -i3 file3.fna
 [-i4 file4.fna]
 [-i5 file5.fna]
 [-i6 file6.fna]
 [-i7 file7.fna]
 [-i8 file8.fna]
 [-i9 file9.fna]
 [-i10 file10.fna]
 [-i11 file11.fna]
 [-i12 file12.fna]
 -orth all_orthomcl.out
 -o orthoset_dir
 -sets 5
 [-funccats1 funccats/pabyss/pabyss_gene_oid_funccat.txt]

=head1 DESCRIPTION

Open a list of 1-1(-1-1-1-1-1-1-1-1-1-1) orthologies and the
multifasta nucleotide files with the ">locus_id_" description, and
create a directory with each set of fastas for each entry in the
orthologies list.

The columns in the orth file and the infiles have to be in the same
order:

ORTHOMCL42(3 genes,3 taxa):      21397376(_NC_003995) 30264936(_NC_003997) 47530432(_NC_007530)

-sets 5   option to indicate which kind of entries are parsed: sets 5
          for a group of 6 genomes will retrieve sets of type
          1-1-1-1-1.

[-altinput] alternative input file for cases when a gene name is a 
            tag in the orth or viceversa (only for 1 file for the moment)

[-funccats1] a comma separated list of (gene_ids,funccats) for input1
             to further subdivide by functional category

=head1 AUTHOR - Albert Vilella

Email 

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut


# Let the code begin...

use strict;
use Bio::SeqIO;
use Getopt::Long;

my ($input1,$input2,$input3,$input4,$input5,
    $input6,$input7,$input8,$input9,$input10,
    $input11,$input12,
    $genes,$taxa,$altinput,$orth,$outdir,$funccats1,
    $spsorder) = undef;
my ($logfile) = undef;
my $merged = undef;

GetOptions(
	   'i1:s' => \$input1,
	   'i2:s' => \$input2,
	   'i3:s' => \$input3,
	   'i4:s' => \$input4,
	   'i5:s' => \$input5,
	   'i6:s' => \$input6,
	   'i7:s' => \$input7,
	   'i8:s' => \$input8,
	   'i9:s' => \$input9,
	   'i10:s' => \$input10,
	   'i11:s' => \$input11,
	   'i12:s' => \$input12,
	   'alt|alti|altinput:s' => \$altinput,
	   'orth:s'  => \$orth,
           'g|genes:s' => \$genes,
           'm|merged' => \$merged,
           't|taxa:s' => \$taxa,
	   'o|outdir:s'=> \$outdir,
	   'f|fun|funccats1:s'=> \$funccats1,
	   'sortsps|spsorder:s'=> \$spsorder,
          );

my $order_by = sort_order("$spsorder");

if (!$outdir) {
    $outdir = "./";
}

unless (-d "$outdir") {
    mkdir "$outdir" or die "couldnt create directory: $!";
}

print "Creating orthomlc alignment sets...\n";
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $logfile_time = 
    sprintf ("%04d%02d%02d_%02d%02d%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec);

$logfile = Bio::Root::IO->new
    (
     -file=>Bio::Root::IO->catfile
     (
      ">","$outdir","logfile.create_orthomcl_alignment_sets.$logfile_time.txt"
     ),
    );

printf ("Writing log to %s\n", $logfile->file);
$logfile->_print
    (
     "Started $logfile_time\n"
    );

my $in1 = new Bio::SeqIO(-file => $input1,
                         -format => 'fasta');
my $in2 = new Bio::SeqIO(-file => $input2,
                         -format => 'fasta');

my $altin='';
if ($altinput) {
    $altin = new Bio::SeqIO(-file => $altinput,
                            -format => 'fasta');
}
my $in3;
my $in4;
my $in5;
my $in6;
my $in7;
my $in8;
my $in9;
my $in10;
my $in11;
my $in12;
if ($input3) {
    $in3 = new Bio::SeqIO(-file => $input3,
                          -format => 'fasta');
    if ($input4) {
        $in4 = new Bio::SeqIO(-file => $input4,
                              -format => 'fasta');
        if ($input5) {
            $in5 = new Bio::SeqIO(-file => $input5,
                                  -format => 'fasta');
            if ($input6) {
                $in6 = new Bio::SeqIO(-file => $input6,
                                      -format => 'fasta');
                if ($input7) {
                    $in7 = new Bio::SeqIO(-file => $input7,
                                          -format => 'fasta');
                    if ($input8) {
                        $in8 = new Bio::SeqIO(-file => $input8,
                                              -format => 'fasta');
                        if ($input9) {
                            $in9 = new Bio::SeqIO(-file => $input9,
                                                  -format => 'fasta');
                            if ($input10) {
                                $in10 = new Bio::SeqIO(-file => $input10,
                                                       -format => 'fasta');
                                if ($input11) {
                                    $in11 = new Bio::SeqIO(-file => $input11,
                                                           -format => 'fasta');
                                    if ($input12) {
                                        $in12 = new Bio::SeqIO(-file => $input12,
                                                               -format => 'fasta');
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

my $count = 0;
my $orthology;
my %hash_of_orthologies;

my $regexp;
my $pattern;

$pattern = 'ORTHOMCL(\S+)\(' . "$genes" . '\ genes\,' . "$taxa" . '\ taxa\)\:';

for ($count=1;$count<=$genes;$count++) {
    $pattern .= '\s+(\S+)\((\S+)\)';
}

$regexp = qr/^${pattern}/i;

open (ORTH,"$orth") or die "cannot open orth file: $!";
print "Loading orthologies...\n";
while (<ORTH>) {
    if ($_ =~ /$regexp/) {
        $orthology = sprintf("%05d", $1);
    $hash_of_orthologies{$orthology}{$3}  = $2;
    $hash_of_orthologies{$orthology}{$5}  = $4;
    eval {$hash_of_orthologies{$orthology}{$7}  =  $6;};
    eval {$hash_of_orthologies{$orthology}{$9}  =  $8;};
    eval {$hash_of_orthologies{$orthology}{$11} = $10;};
    eval {$hash_of_orthologies{$orthology}{$13} = $12;};
    eval {$hash_of_orthologies{$orthology}{$15} = $14;};
    eval {$hash_of_orthologies{$orthology}{$17} = $16;};
    eval {$hash_of_orthologies{$orthology}{$19} = $18;};
    eval {$hash_of_orthologies{$orthology}{$21} = $20;};
    eval {$hash_of_orthologies{$orthology}{$23} = $22;};
    eval {$hash_of_orthologies{$orthology}{$25} = $24;};
    $orthology = "";
    $count++;
    }
}

close ORTH;

print "Loaded $count orthologous entries...\n";

my %seqs = undef;
my ($display_id);
print "Loading file1...\n";
my $num_processed = 0;
while ( my $seq = $in1->next_seq ) {
    $logfile->_print("file $input1\n");
    $logfile->_print("no_$num_processed $display_id\n");
    $num_processed++;
    $display_id = $seq->display_id;
    $seqs{$display_id} = $seq;
}

$num_processed=0;
print "Loading file2...\n";
while ( my $seq = $in2->next_seq ) {
    $logfile->_print("file $input2\n");
    $logfile->_print("no_$num_processed $display_id\n");
    $num_processed++;
    $display_id = $seq->display_id;
    $seqs{$display_id} = $seq;
}

$num_processed=0;
print "Loading file3...\n";
if ($input3) {
    while ( my $seq = $in3->next_seq ) {
        $logfile->_print("file $input3 \n");
        $logfile->_print("no_$num_processed $display_id\n");
        $num_processed++;
        $display_id = $seq->display_id;
        $seqs{$display_id} = $seq;
    }
}

$num_processed=0;
if ($altinput) {
    print "Loading altseq...\n";
    while ( my $seq = $altin->next_seq ) {
        $logfile->_print("file $altinput\n");
        $logfile->_print("no_$num_processed $display_id\n");
        $num_processed++;
        $display_id = $seq->display_id;
        $seqs{$display_id} = $seq;
    }
}

$num_processed=0;
if ($input4) {
    print "Loading file4...\n";
    while ( my $seq = $in4->next_seq ) {
        $logfile->_print("file $input4 \n");
        $logfile->_print("no_$num_processed $display_id\n");
        $num_processed++;
        $display_id = $seq->display_id;
        $seqs{$display_id} = $seq;
    }
}

$num_processed=0;
if ($input5) {
    print "Loading file5...\n";
    while ( my $seq = $in5->next_seq ) {
        $logfile->_print("file $input5 \n");
        $logfile->_print("no_$num_processed $display_id\n");
        $num_processed++;
        $display_id = $seq->display_id;
        $seqs{$display_id} = $seq;
    }
}

$num_processed=0;
if ($input6) {
    print "Loading file6...\n";
    while ( my $seq = $in6->next_seq ) {
        $logfile->_print("file $input6 \n");
        $logfile->_print("no_$num_processed $display_id\n");
        $num_processed++;
        $display_id = $seq->display_id;
        $seqs{$display_id} = $seq;
    }
}

$num_processed=0;
if ($input7) {
    print "Loading file7...\n";
    while ( my $seq = $in7->next_seq ) {
        $logfile->_print("file $input7 \n");
        $logfile->_print("no_$num_processed $display_id\n");
        $num_processed++;
        $display_id = $seq->display_id;
        $seqs{$display_id} = $seq;
    }
}

$num_processed=0;
if ($input8) {
    print "Loading file8...\n";
    while ( my $seq = $in8->next_seq ) {
        $logfile->_print("file $input8 \n");
        $logfile->_print("no_$num_processed $display_id\n");
        $num_processed++;
        $display_id = $seq->display_id;
        $seqs{$display_id} = $seq;
    }
}

$num_processed=0;
if ($input9) {
    print "Loading file9...\n";
    while ( my $seq = $in9->next_seq ) {
        $logfile->_print("file $input9 \n");
        $logfile->_print("no_$num_processed $display_id\n");
        $num_processed++;
        $display_id = $seq->display_id;
        $seqs{$display_id} = $seq;
    }
}

$num_processed=0;
if ($input10) {
    print "Loading file10...\n";
    while ( my $seq = $in10->next_seq ) {
        $logfile->_print("file $input10 \n");
        $logfile->_print("no_$num_processed $display_id\n");
        $num_processed++;
        $display_id = $seq->display_id;
        $seqs{$display_id} = $seq;
    }
}

$num_processed=0;
if ($input11) {
    print "Loading file11...\n";
    while ( my $seq = $in11->next_seq ) {
        $logfile->_print("file $input11\n");
        $logfile->_print("no_$num_processed $display_id\n");
        $num_processed++;
        $display_id = $seq->display_id;
        $seqs{$display_id} = $seq;
    }
}

$num_processed=0;
if ($input12) {
    print "Loading file12...\n";
    while ( my $seq = $in12->next_seq ) {
        $logfile->_print("file $input12\n");
        $logfile->_print("no_$num_processed $display_id\n");
        $num_processed++;
        $display_id = $seq->display_id;
        $seqs{$display_id} = $seq;
    }
}

$num_processed=0;

my %funccats1 = undef;
$count = 0;
if ($funccats1) {
    open (FUNCCATS1,"$funccats1") or die "cannot open orth file: $!";
    print "Loading funccats1...\n";
    while (<FUNCCATS1>) {
        if ($_ =~ /(\S+)\,(\S+)/g) {
            $funccats1{$1} = $2;
            $count++;
        } else {
            warn "error parsing funccats1: $_\n";
        }
    }
}
print "Loaded $count functional category entries...\n";

my $out;
my $outseq;
my $geneid = undef;
$count = 0;
print "Creating sets...\n";
foreach my $orthonumber (sort keys %hash_of_orthologies) {
    my $spnum = 1; my $firstgene = "";
    my @list = keys %{ $hash_of_orthologies{$orthonumber} };
    @list = sort { $order_by->{$a} <=> $order_by->{$b} } @list;
    foreach my $species (@list) {
#     foreach my $species (sort keys %{ $hash_of_orthologies{$orthonumber} }) {
        $geneid = $hash_of_orthologies{$orthonumber}{$species};
        if (length($geneid) == 0) {
            next;
        } # get rid of null evals
        if ($spnum == 1) {
            $firstgene = $geneid;
            my $subdir = undef;
            if ($funccats1) {
                my $cat = $funccats1{$firstgene} || "w";
                unless (-d "$outdir/$cat") {
                    mkdir ("$outdir/$cat");
                }
                $subdir = "$outdir/$cat/$firstgene";
            } else {
                $subdir = "$outdir/$firstgene";
            }
            unless (-d "$subdir") {
                mkdir ("$subdir");
            }
            $logfile->_print
                (
                 "$subdir\n"
                );
            $out = new Bio::SeqIO(-file => ">$subdir/$geneid.fasta");
        }
        if (($seqs{$geneid})) {
            $outseq = $seqs{$geneid};
            my $species_mark_primary_id = "sp" . sprintf("%02d", $spnum);
            $species_mark_primary_id .= "_" . $species . "_" . $outseq->primary_id;
            $outseq->display_id("$species_mark_primary_id");
            $out->write_seq($outseq);
        } else {
            warn "have not found $geneid in $species for orthology $orthonumber\n";
        }
        $spnum++;
    }
}

$logfile->_print("End\n");
$logfile->close();

print "\n";
print "Done.\n";

1;

########################################

sub sort_order {
    my @names = split /\:/,shift;
    my %order_by = map { $names[$_] => $_ } 0 .. $#names;
    if ( keys %order_by ne @names ) {
        warn "duplicate names found";
    }
    return \%order_by;
}


1;
