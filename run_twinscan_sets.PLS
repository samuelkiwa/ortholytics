#!/usr/local/bin/perl
# run_twinscan_sets.PLS
#
# Cared for by Filipe G. Vieira <>
#
# Copyright Filipe G. Vieira
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

create_twinscan_sets.PLS - DESCRIPTION

=head1 SYNOPSIS

perl run_twinscan_sets.PLS \
-refsp dmel \
-targetsp dsec \
-exedir /oath/to/where/twinscan/is/installed
-seqdir /path/to/both/melano/and/sechellia/GFF/and/fasta/files \
-outdir /path/to/where/each/set/subdir/should/be/created \
-twin (if set, does not create de novo FASTA files)

=head1 DESCRIPTION

It script will take a set of seqeunces in FASTA of two species
(ref and target) and perform an all-against-all prediction
of genes using Twinscan.

If -twin option is used, the script will assume that the sequences already
exist and will, therefore, just run twinscan. Note that the sequences
must be in the outdir/seq/ directory.

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
use File::Copy;
use File::Basename;
use Bio::SeqIO;
use Bio::Root::IO;

use constant FLANK_SIZE => 200;

my ($seq_dir, $out_dir, $out_name, $out_seqs, $out_fasta);
my ($seq, $ref_sp, $target_sp, $exe_dir, $rm_dir,
    $in_ref_fasta, $in_target_fasta, $twin_flag,$onlycreateseqsets,$onlyrepeatmasker);
my (@files, @ref_seqs, @target_seqs);

#Get the command-line options
GetOptions
    (
     'r|refsp:s'  => \$ref_sp,       #reference species
     't|targetsp:s' => \$target_sp,  #target species
     'e|exedir:s' => \$exe_dir,      #directory where Twinscan is installed
     's|seqdir:s' => \$seq_dir,      #directory to search for files
     'o|outdir:s' => \$out_dir,      #directory to output
     'w|twin' => \$twin_flag,         #just run twinsscan
     'c|createsets|only' => \$onlycreateseqsets,         #just create seq sets
     'onlyrepeatmasker' => \$onlyrepeatmasker,         #just call repeatmasker and leave
);

$out_seqs = Bio::Root::IO->catfile($out_dir,"seqs");

goto TWINSCAN if($twin_flag);
mkpath($out_seqs);

@files="";
find( sub{ m/^$ref_sp-all-chromosome-r[0-9].+\.fasta$/i and
  unshift(@files,$File::Find::name);}, $seq_dir, );
if ($files[0] eq ""){print("$ref_sp FASTA files not found.\n");exit(0)}
$in_ref_fasta = new Bio::SeqIO(-file => $files[0], -format => "Fasta");
print STDERR "Loading ref_sp fasta files...\n";
while($seq = $in_ref_fasta->next_seq)
{
    $out_name = Bio::Root::IO->catfile($out_seqs, $ref_sp."_".$seq->display_id.".fasta");
    $out_fasta = new Bio::SeqIO(-file => ">".$out_name, -format => "Fasta");
    $out_fasta->write_seq($seq);
    $out_fasta->close;
}
$in_ref_fasta->close;


@files="";
find( sub{ m/^$target_sp\.scaffolds\.fasta$/i and
  unshift(@files,$File::Find::name);}, $seq_dir, );
if ($files[0] eq ""){print("$target_sp FASTA files not found.\n");exit(0)}
$in_target_fasta = new Bio::SeqIO(-file => $files[0], -format => "Fasta");
print STDERR "Loading target_sp fasta files...\n";
while($seq = $in_target_fasta->next_seq)
{
    $out_name = Bio::Root::IO->catfile($out_seqs, $target_sp."_".$seq->display_id.".fasta");
    $out_fasta = new Bio::SeqIO(-file => ">".$out_name, -format => "Fasta");
    $out_fasta->write_seq($seq);
    $out_fasta->close;
}
$in_target_fasta->close;

print STDERR "Created seq sets...\n";
exit if ($onlycreateseqsets);

TWINSCAN:

print STDERR "All seqs set. Starting twinscan loop...\n";

@files="";
find( sub{ m/\.fasta$/i and
  unshift(@files,$File::Find::name);}, $out_seqs, );

@ref_seqs = grep(m/$ref_sp/i, @files);
@target_seqs = grep(m/$target_sp/i, @files);

my @sort;
@sort = sort {$a <=> $b} @ref_seqs;
@ref_seqs = @sort;

@sort = sort {$a <=> $b} @target_seqs;
@target_seqs = @sort;


foreach my $ref_file (@ref_seqs)
{
    my ($ref_chr_name, $target_chr_name, $ref_dir, $out_file);
    
    $ref_chr_name = basename($ref_file);
    $ref_chr_name =~ s/\.fasta//g;
    $ref_dir = Bio::Root::IO->catfile($out_dir,$ref_chr_name);
    mkpath($ref_dir);
    
    foreach my $target_file (@target_seqs)
    {
        $target_chr_name = basename($target_file);
        $out_file = Bio::Root::IO->catfile($out_dir,$target_chr_name);
	print($ref_chr_name." <=> ".$target_chr_name."\n");

        ##Run twinscan
        my $exe = Bio::Root::IO->catfile($exe_dir, "bin", "Twinscan.pl");
        my $zoe = Bio::Root::IO->catfile($exe_dir, "parameters", "iscan.dm2.albert_vilella.zhmm");
        my $call =
            "perl $exe "
            ."-r $zoe "
            ."-d $out_dir ";
        $call .= "-onlyrepeatmasker " if ($onlyrepeatmasker);
        $call .= 
            "$target_file "
            ."$ref_file "
            ."1>>$out_dir/run_twin.log "
            ."2>>$out_dir/run_twin.err";
        print STDERR "Calling this now:\n";
        print STDERR "$call\n";
        system($call);

        move($out_file.".blast", $ref_dir);
        move($out_file.".cat", $ref_dir);
        move($out_file.".conseq", $ref_dir);
        move($out_file.".gff", $ref_dir);
        move($out_file.".out", $ref_dir);
        move($out_file.".tbl", $ref_dir);
        move($out_file.".twinscan", $ref_dir);
        move($target_file.".twinscan.progress", $ref_dir);
    }
}

1;

