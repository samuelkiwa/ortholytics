#!/usr/local/bin/perl
# create_projector_sets.PLS
#
# Cared for by Filipe G. Vieira <>
#
# Copyright Filipe G. Vieira
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

create_projector_sets.PLS - DESCRIPTION 

=head1 SYNOPSIS

perl create_projector_sets.PLS \
-refsp dmel \
-targetsp dsec \
-seqdir /path/to/both/melano/and/sechellia/GFF/and/fasta/files \
-outdir /path/to/where/each/set/subdir/should/be/created \

=head1 DESCRIPTION

This script will take a template target (or known) gff file, for
example, melano, and filter in each entry that relates to CDS (genes
and exons).

It will also take a query (genes not known) gff file from a tblastn
run for the assembled species (for example, sechellia) and use the
start and end plus 200bp up+downstream as query sequences.

It finally take the fasta seqs in the fastadir option and, jointly
with the gff information, create each subdir for each Projector
execution inside outdir.

After this script, all should be ready for running Projector in each
subdir.

=head1 AUTHOR - Filipe G. Vieira

Email 

Describe contact details here

=head1 CONTRIBUTORS

Albert Vilella

=cut

# Let the code begin...

use strict;
use Getopt::Long;
use File::Path;
use File::Find;
use POSIX;
use Bio::SeqIO;
use Bio::Root::IO;
use Bio::Tools::GFF;

use constant FLANK_SIZE => 200;

my ($ref_sp, $target_sp,
    $in_ref_fasta, $in_ref_gff,
    $in_target_gff, $in_target_fasta,
    $dir, $outdir);
my ($seq, @files, 
    %ref_chr, %target_chr,
    %ref_genes, %target_genes);
my ($rel_pos, $exon_n);

#Get the command-line options
GetOptions(
           'r|refsp:s'  => \$ref_sp,       #reference species
           't|targetsp:s' => \$target_sp,  #target species
           's|seqdir:s' => \$dir,          #directory to search for files
           'o|outdir:s' => \$outdir        #directory to output
	    );

@files="";
find( sub{ m/^$ref_sp-all-chromosome-r[0-9].+\.fasta$/i and
    unshift(@files,$File::Find::name);}, $dir, );
if ($files[0] eq ""){print("$ref_sp FASTA files not found.\n");exit(0)}
$in_ref_fasta = new Bio::SeqIO(-file => $files[0], -format => "largefasta");
#Load ref_sp Fasta file
while($seq = $in_ref_fasta->next_seq)
{$ref_chr{$seq->display_id} = $seq;}

@files="";
find( sub{ m/^$target_sp\.scaffolds\.fasta$/i and
    unshift(@files,$File::Find::name);}, $dir, );
if ($files[0] eq ""){print("$target_sp FASTA files not found.\n");exit(0)}
$in_target_fasta = new Bio::SeqIO(-file => $files[0], -format => "Fasta");
#Load target_sp Fasta file
while($seq = $in_target_fasta->next_seq)
{$target_chr{$seq->display_id} = $seq;}

@files="";
find( sub{ m/^$ref_sp-.+-r[0-9].+\.gff$/i and
    unshift(@files,$File::Find::name);}, $dir, );
if ($files[0] eq ""){print("$ref_sp GFF files not found.\n");exit(0)}
foreach my $file (@files)
{
    $in_ref_gff = new Bio::Tools::GFF
        (
         -file => $file, 
         -gff_version => 3,
        );
    while ($seq = $in_ref_gff->next_feature)
    {
        my ($id, @tags);
        if($seq->primary_tag eq "exon")
        {
            @tags = $seq->get_tag_values("ID");
            @tags = split(m/:/, $tags[0]);
            $id = $tags[1];
            @tags = $seq->get_tag_values("Parent");
            foreach my $parent (@tags)
            {
                $ref_genes{$parent}{"exon"}[$id]{"start"} = $seq->start;
                $ref_genes{$parent}{"exon"}[$id]{"end"} = $seq->end;
                $ref_genes{$parent}{"exon"}[$id]{"feature"} = $seq;
            }
        }
        elsif($seq->primary_tag eq "CDS")
        {
            @tags = $seq->get_tag_values("Parent");
            $ref_genes{$tags[0]}{"CDS"}[0]{"chr"} = $seq->seq_id;
            $ref_genes{$tags[0]}{"CDS"}[0]{"start"} = $seq->start;
            $ref_genes{$tags[0]}{"CDS"}[0]{"end"} = $seq->end;
            $ref_genes{$tags[0]}{"CDS"}[0]{"strand"} =
   		 (($seq->strand == 1)?("+"):("-"));
        }
    }
    $in_ref_gff->close;
}

#Clear empty exon postitions and sort exons
foreach $seq (keys(%ref_genes))
{
    my @sort;
    @{$ref_genes{$seq}{"exon"}} = 
    	grep(defined, @{$ref_genes{$seq}{"exon"}});

    if($ref_genes{$seq}{"CDS"}[0]{"strand"} eq "+")
    {@sort = sort {$a->{"start"} <=> $b->{"start"}} @{$ref_genes{$seq}{"exon"}};}
    else
    {@sort = reverse sort {$a->{"start"} <=> $b->{"start"}}
          @{$ref_genes{$seq}{"exon"}};}

    @{$ref_genes{$seq}{"exon"}} = @sort;
    unshift(@{$ref_genes{$seq}{"exon"}}, undef);
}

foreach $seq (sort keys(%ref_genes))
{
    for($exon_n=1; $ref_genes{$seq}{"exon"}[$exon_n] ne undef; $exon_n++)
    {
        if ($ref_genes{$seq}{"CDS"}[0]{"start"} <= $ref_genes{$seq}{"exon"}[$exon_n]{"end"} &&
            $ref_genes{$seq}{"CDS"}[0]{"end"} >= $ref_genes{$seq}{"exon"}[$exon_n]{"start"})
        {
        	  $ref_genes{$seq}{"CDS"}[$exon_n]{"start"} = $ref_genes{$seq}{"exon"}[$exon_n]{"start"};
        	  $ref_genes{$seq}{"CDS"}[$exon_n]{"end"} = $ref_genes{$seq}{"exon"}[$exon_n]{"end"};

        	  $ref_genes{$seq}{"CDS"}[$exon_n]{"start"} = $ref_genes{$seq}{"CDS"}[0]{"start"}
                if ($ref_genes{$seq}{"CDS"}[0]{"start"} > $ref_genes{$seq}{"exon"}[$exon_n]{"start"});
        	  $ref_genes{$seq}{"CDS"}[$exon_n]{"end"} = $ref_genes{$seq}{"CDS"}[0]{"end"}
                if ($ref_genes{$seq}{"CDS"}[0]{"end"} < $ref_genes{$seq}{"exon"}[$exon_n]{"end"});
            $ref_genes{$seq}{"CDS"}[$exon_n]{"exon_n"} = $exon_n;
        }
    }
}

foreach $seq (keys(%ref_genes))
{@{$ref_genes{$seq}{"CDS"}} = grep(defined, @{$ref_genes{$seq}{"CDS"}});}

@files="";
find( sub{ m/^$target_sp-prot9-strong2\.gff$/i and
    unshift(@files,$File::Find::name);}, $dir, );
if ($files[0] eq ""){print("$target_sp GFF files not found.\n");exit(0)}
foreach my $file (@files)
{
    $in_target_gff = new Bio::Tools::GFF
        (
         -file => $file,
         -gff_version => 3,
        );

    while ($seq = $in_target_gff->next_feature)
    {
        my ($id, $start, $end, $len, @tags);

        if($seq->source_tag eq "tblastn:modDM" && $seq->primary_tag eq "match")
        {
            @tags = $seq->get_tag_values("ID");
            $id = $tags[0];
            $id =~ s/-P/-R/g;

            next unless(defined($ref_genes{$id}));

            my $output_dir = Bio::Root::IO->catfile($outdir,"proj_".$id);
            mkpath($output_dir);

            open(OUT_FASTA, ">".Bio::Root::IO->catfile($output_dir,$target_sp.".".$id.".fasta"));
            #get ref_sp sequence
            $start = $ref_genes{$id}{"CDS"}[0]{"start"} - FLANK_SIZE;
            $end = $ref_genes{$id}{"CDS"}[0]{"end"} + FLANK_SIZE;
            if($start < 0){$start = 1;}
            if($end > $ref_chr{$ref_genes{$id}{"CDS"}[0]{"chr"}}->length)
	    {$end = $ref_chr{$ref_genes{$id}{"CDS"}[0]{"chr"}}->length;}
           $rel_pos = $start - 1; #pos to update GTF

            #Print reference sequence info
            print($id." => ".$ref_genes{$id}{"CDS"}[0]{"chr"}
		  ."/".$ref_genes{$id}{"CDS"}[0]{"strand"}
		  ."/".$start."-".$end." <=> ");

            my $subseq = $ref_chr{$ref_genes{$id}{"CDS"}[0]{"chr"}}->subseq($start,$end);
            $len = length($subseq);
            print(OUT_FASTA ">".$id." 1-".$len." "
      .(($ref_genes{$id}{"CDS"}[0]{"strand"} eq "+")?("forward"):("reverse"))."\n");
            print(OUT_FASTA $subseq."\n");

            #get target_sp sequence
            my $chr = $seq->seq_id;
            $chr =~ s/^scaffold_/super_/g;
            $start = floor($seq->start + ($seq->end - $seq->start)/2 - $len/2 - FLANK_SIZE);
            $end = ceil($start + $len + 2*FLANK_SIZE);
            if($start < 0){$start = 1;}
            if($end > $target_chr{$chr}->length){$end = $target_chr{$chr}->length;}

            #Print target sequence info
            print($chr."/".(($seq->strand eq "1")?("+"):("-"))
		  ."/".$start."-".$end."\n");

            my $subseq = $target_chr{$chr}->subseq($start,$end);
            print(OUT_FASTA ">".$target_sp."$id 1-".length($subseq)." "
		  .(($seq->strand eq "1")?("forward"):("reverse"))."\n");
            print(OUT_FASTA $subseq."\n");
            close(OUT_FASTA);

            open(OUT_GTF, ">".Bio::Root::IO->catfile($output_dir,
			    			     $ref_sp.".".$id.".gtf"));

            my @gtf = &gff2gtf($rel_pos, $id, \%ref_genes);
            print(OUT_GTF join("\n",@gtf));
            close(OUT_GTF);
        }
    }
}

########################################
## FUNCTIONS
########################################

sub gff2gtf()
{
    my ($pos, $id, $genes) = @_;

    my ($gene_name, $transc_name, $exon_n, $cds_exon,
        $start, $end, $phase, $strand, @gtf, $prev_phase, %genes);

    my ($start_codon_pos, $stop_codon_pos) = (3,3);
    my ($cds_flag, $cds_exon) = (0,1);

    next unless ($id =~ m/^CG|^GA/);
    $gene_name = substr($id,0,-3);
    $transc_name = $id;
    $strand = $genes->{$id}{"CDS"}[0]{"strand"};

    for($exon_n=1, $cds_exon=1; $genes->{$id}{"exon"}[$exon_n] ne undef; $exon_n++)
    {
        ########################################################
        #Add GTF features
        foreach my $type ("exon", "start_codon", "CDS", "stop_codon")
        {
            if($type eq "exon")
            {
                $start = $genes->{$id}{"exon"}[$exon_n]{"start"};
                $end = $genes->{$id}{"exon"}[$exon_n]{"end"};
                $phase = ".";
            }

            elsif($type eq "start_codon" &&
                  $start_codon_pos &&
                  $exon_n >= $genes->{$id}{"CDS"}[1]{"exon_n"})
            {
                if($strand eq "+")
                {
                    $start = $genes->{$id}{"CDS"}[$cds_exon]{"start"};
                    $end = $start + ($start_codon_pos - 1);
                    $end = $genes->{$id}{"CDS"}[$cds_exon]{"end"}
                        if ($end > $genes->{$id}{"CDS"}[$cds_exon]{"end"});
                }
                else
                {
                    $start = $genes->{$id}{"CDS"}[$cds_exon]{"end"} - ($start_codon_pos - 1);
                    $end = $genes->{$id}{"CDS"}[$cds_exon]{"end"};
                    $start = $genes->{$id}{"CDS"}[$cds_exon]{"start"}
                        if ($start < $genes->{$id}{"CDS"}[$cds_exon]{"start"});
                }
                $phase = "0";
                $start_codon_pos -= $end - $start + 1;
            }

            elsif($type eq "stop_codon" &&
                  $stop_codon_pos &&
                  $exon_n >= $genes->{$id}{"CDS"}[-1]{"exon_n"})
            {
                if($strand eq "+")
                {
                    if($genes->{$id}{"CDS"}[$cds_exon]{"end"} ne "")
                    {$start = $genes->{$id}{"CDS"}[$cds_exon]{"end"} + 1;}
                    else
                    {$start = $genes->{$id}{"exon"}[$exon_n]{"start"};}

                    $end = $start + $stop_codon_pos - 1;
                    $end = $genes->{$id}{"exon"}[$exon_n]{"end"}
                        if ($end > $genes->{$id}{"exon"}[$exon_n]{"end"});
                    next if ($start > $end);
                }
                else
                {
                    if($genes->{$id}{"CDS"}[$cds_exon]{"start"} ne "")
                    {$start = $genes->{$id}{"CDS"}[$cds_exon]{"start"} - $stop_codon_pos;}
                    else
                    {$start = $genes->{$id}{"exon"}[$exon_n]{"end"} - ($stop_codon_pos - 1);}

                    $end = $start + ($stop_codon_pos - 1);

                    $start = $genes->{$id}{"exon"}[$exon_n]{"start"}
                        if ($start < $genes->{$id}{"exon"}[$exon_n]{"start"});
                    next if ($start > $end);
                }
                $phase = "0";
                $stop_codon_pos -= $end - $start + 1;
            }

            elsif($type eq "CDS" &&
                  $exon_n == $genes->{$id}{"CDS"}[$cds_exon]{"exon_n"})
            {
                $cds_flag = 1;
                $start = $genes->{$id}{"CDS"}[$cds_exon]{"start"};
                $end = $genes->{$id}{"CDS"}[$cds_exon]{"end"};

                $phase = (($prev_phase == 0)?(0):(3 - $prev_phase));
                $prev_phase = ($end-$start+1-$phase)%3;
            }
            else{next;}

            push(@gtf, sprintf("%s\tFlyBase\t%s\t%s\t%s\t.\t%s\t%s\tgene_id \"%s\"; transcript_id \"%s\"; exon_number \"%s\";",
                   $id,
                   $type,
                   $start - $pos,
                   $end - $pos,
                   $strand,
                   $phase,
                   $gene_name,
                   $transc_name,
                   $exon_n));
        }
        if($cds_flag){$cds_exon++;}
    }

    return (@gtf);
}

1;
