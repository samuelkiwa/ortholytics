#!/usr/local/bin/perl
# gff2gtf.PLS
#
# Cared for by Filipe G. Vieira <>
#
# Copyright Filipe G. Vieira
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

gff2gtf.PLS - DESCRIPTION

=head1 SYNOPSIS

perl create_projector_sets.PLS \
-file GFF file to convert to GTF \

=head1 DESCRIPTION

This script will read a GFF file and convert it to the GTF format.

=head1 AUTHOR - Filipe G. Vieira

Email

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut

# Let the code begin...

use strict;
use Getopt::Long;
use Bio::Tools::GFF;

my ($file_name, $gff_file, $id, $seq, @tags, %genes);
my ($gene_name, $transc_name, $exon_n,
    $cds_exon, $start, $end, $phase, $strand, $prev_phase);

#Get the command-line options
&GetOptions(
      "f|file:s" => \$file_name         #directory to search for files
);


$gff_file = new Bio::Tools::GFF(-file => $file_name, -gff_version => 3);

while ($seq = $gff_file->next_feature)
{
    if($seq->primary_tag eq "exon")
    {
        @tags = $seq->get_tag_values("ID");
        @tags = split(m/:/, $tags[0]);
        $id = $tags[1];
        @tags = $seq->get_tag_values("Parent");
        foreach my $parent (@tags)
        {
            $genes{$parent}{"exon"}[$id]{"start"} = $seq->start;
            $genes{$parent}{"exon"}[$id]{"end"} = $seq->end;
            $genes{$parent}{"exon"}[$id]{"feature"} = $seq;
        }
    }
    elsif($seq->primary_tag eq "CDS")
    {
        @tags = $seq->get_tag_values("Parent");
        $genes{$tags[0]}{"CDS"}[0]{"chr"} = $seq->seq_id;
        $genes{$tags[0]}{"CDS"}[0]{"start"} = $seq->start;
        $genes{$tags[0]}{"CDS"}[0]{"end"} = $seq->end;
        $genes{$tags[0]}{"CDS"}[0]{"strand"} = (($seq->strand == 1)?("+"):("-"));
    }
}

$gff_file->close;

#Clear empty exon postitions and sort exons
foreach $seq (keys(%genes))
{
    my @sort;
    #grep defined and sort exons
    @{$genes{$seq}{"exon"}} = grep(defined, @{$genes{$seq}{"exon"}});
    @sort = sort {$a->{"start"} <=> $b->{"start"}} @{$genes{$seq}{"exon"}};
    if($genes{$seq}{"CDS"}[0]{"strand"} eq "-"){@sort = reverse (@sort);}
    @{$genes{$seq}{"exon"}} = @sort;
    unshift(@{$genes{$seq}{"exon"}}, undef);
}


foreach $seq (sort keys(%genes))
{
    for($exon_n=1; $genes{$seq}{"exon"}[$exon_n] ne undef; $exon_n++)
    {
        if ($genes{$seq}{"CDS"}[0]{"start"} <= $genes{$seq}{"exon"}[$exon_n]{"end"} &&
            $genes{$seq}{"CDS"}[0]{"end"} >= $genes{$seq}{"exon"}[$exon_n]{"start"})
        {
        	  $genes{$seq}{"CDS"}[$exon_n]{"start"} = $genes{$seq}{"exon"}[$exon_n]{"start"};
        	  $genes{$seq}{"CDS"}[$exon_n]{"end"} = $genes{$seq}{"exon"}[$exon_n]{"end"};
        
        	  $genes{$seq}{"CDS"}[$exon_n]{"start"} = $genes{$seq}{"CDS"}[0]{"start"}
                  if ($genes{$seq}{"CDS"}[0]{"start"} > $genes{$seq}{"exon"}[$exon_n]{"start"});
        	  $genes{$seq}{"CDS"}[$exon_n]{"end"} = $genes{$seq}{"CDS"}[0]{"end"}
                  if ($genes{$seq}{"CDS"}[0]{"end"} < $genes{$seq}{"exon"}[$exon_n]{"end"});
            $genes{$seq}{"CDS"}[$exon_n]{"exon_n"} = $exon_n;
        }
    }
}

foreach $seq (keys(%genes))
{@{$genes{$seq}{"CDS"}} = grep(defined, @{$genes{$seq}{"CDS"}});}

$file_name =~ s/\.gff$/\.gtf/g;
open(GTF, ">".$file_name) or die ("Error: could not create GTF file.");
foreach $seq (sort keys(%genes))
{
    my ($start_codon_pos, $stop_codon_pos) = (3,3);
    my ($cds_flag, $cds_exon) = (0,1);
    
    next unless ($seq =~ m/^CG|^GA/);
    $gene_name = substr($seq,0,-3);
    $transc_name = $seq;
    $strand = $genes{$seq}{"CDS"}[0]{"strand"};

    for($exon_n=1, $cds_exon=1; $genes{$seq}{"exon"}[$exon_n] ne undef; $exon_n++)
    {
        ########################################################
        #Add GTF features
        foreach my $type ("exon", "start_codon", "CDS", "stop_codon")
        {
            if($type eq "exon")
            {
                $start = $genes{$seq}{"exon"}[$exon_n]{"start"};
                $end = $genes{$seq}{"exon"}[$exon_n]{"end"};
                $phase = ".";
            }

            elsif($type eq "start_codon" &&
                  $start_codon_pos &&
                  $exon_n >= $genes{$seq}{"CDS"}[1]{"exon_n"})
            {
                if($strand eq "+")
                {
                    $start = $genes{$seq}{"CDS"}[$cds_exon]{"start"};
                    $end = $start + ($start_codon_pos - 1);
                    $end = $genes{$seq}{"CDS"}[$cds_exon]{"end"}
                          if ($end > $genes{$seq}{"CDS"}[$cds_exon]{"end"});
                }
                else
                {
                    $start = $genes{$seq}{"CDS"}[$cds_exon]{"end"} - ($start_codon_pos - 1);
                    $end = $genes{$seq}{"CDS"}[$cds_exon]{"end"};
                    $start = $genes{$seq}{"CDS"}[$cds_exon]{"start"}
                          if ($start < $genes{$seq}{"CDS"}[$cds_exon]{"start"});
                }
                $phase = "0";
                $start_codon_pos -= $end - $start + 1;
            }

            elsif($type eq "stop_codon" &&
                  $stop_codon_pos &&
                  $exon_n >= $genes{$seq}{"CDS"}[-1]{"exon_n"})
            {
                if($strand eq "+")
                {
                    if($genes{$seq}{"CDS"}[$cds_exon]{"end"} ne "")
                    {$start = $genes{$seq}{"CDS"}[$cds_exon]{"end"} + 1;}
                    else
                    {$start = $genes{$seq}{"exon"}[$exon_n]{"start"};}
                    
                    $end = $start + $stop_codon_pos - 1;
                    $end = $genes{$seq}{"exon"}[$exon_n]{"end"}
                          if ($end > $genes{$seq}{"exon"}[$exon_n]{"end"});
                    next if ($start > $end);
                }
                else
                {
                    if($genes{$seq}{"CDS"}[$cds_exon]{"start"} ne "")
                    {$start = $genes{$seq}{"CDS"}[$cds_exon]{"start"} - $stop_codon_pos;}
                    else
                    {$start = $genes{$seq}{"exon"}[$exon_n]{"end"} - ($stop_codon_pos - 1);}

                    $end = $start + ($stop_codon_pos - 1);

                    $start = $genes{$seq}{"exon"}[$exon_n]{"start"}
                          if ($start < $genes{$seq}{"exon"}[$exon_n]{"start"});
                    next if ($start > $end);
                }
                $phase = "0";
                $stop_codon_pos -= $end - $start + 1;
            }

            elsif($type eq "CDS" &&
                  $exon_n == $genes{$seq}{"CDS"}[$cds_exon]{"exon_n"})
            {
                $cds_flag = 1;
                $start = $genes{$seq}{"CDS"}[$cds_exon]{"start"};
                $end = $genes{$seq}{"CDS"}[$cds_exon]{"end"};

                $phase = (($prev_phase == 0)?(0):(3 - $prev_phase));
                $prev_phase = ($end-$start+1-$phase)%3;
            }
            else{next;}

            printf(GTF "%s\tFlyBase\t%s\t%s\t%s\t.\t%s\t%s\tgene_id \"%s\"; transcript_id \"%s\"; exon_number \"%s\";\n",
                   $genes{$seq}{"exon"}[$exon_n]{"feature"}->seq_id,
                   $type,
                   $start,
                   $end,
                   $strand,
                   $phase,
                   $gene_name,
                   $transc_name,
                   $exon_n);
        }
        if($cds_flag){$cds_exon++;}
    }
}
close (GTF);
print("GTF successfully created!\n");

1;


