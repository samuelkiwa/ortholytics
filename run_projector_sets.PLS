#!/usr/local/bin/perl
# run_projector_sets.PLS
#
# Cared for by Filipe G. Vieira <>
#
# Copyright Filipe G. Vieira
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

run_projector_sets.PLS - DESCRIPTION

=head1 SYNOPSIS

perl run_projector_sets.PLS \
-seqdir /path/to/where/each/set/subdir/is \
-exedir /path/to/where/projector/is/isntalled

=head1 DESCRIPTION

This script will run projector on a directory tree previously
created by "run_projector_sets.PLS".

=head1 AUTHOR - Filipe G. Vieira

Email

Describe contact details here

=head1 CONTRIBUTORS

Albert Vilella

=cut


# Let the code begin...

use strict;
use Getopt::Long;
use File::Find;
use File::Copy;
use File::Basename;
use Bio::Root::IO;

my ($seq_dir, $exe_dir, @tmp_file, @exec_files, @files);

#Get the command-line options
GetOptions(
      's|seqdir=s' => \$seq_dir,         #directory to search for files
      'e|exedir=s' => \$exe_dir,         #directory of projector executable
);

find( sub{ m/\.fasta$|\.gtf$/i and unshift(@files,$File::Find::name);}, $seq_dir, );
@files = sort @files;

do{
    @exec_files = splice(@files, 0, 2);

    @tmp_file = grep(m/\.fasta$/g, @exec_files);
    my $data_dir = dirname($tmp_file[0]);
    my $fasta_file = basename($tmp_file[0]);
    
    @tmp_file = grep(m/\.gtf$/g, @exec_files);
    my $gtf_file = basename($tmp_file[0]);

    my $prefix = $fasta_file;
    $prefix =~ s/\.fasta//g;

    my $exe = Bio::Root::IO->catfile($exe_dir,"projector2.0beta");
    my $log = Bio::Root::IO->catfile($data_dir,"run_proj.log");
    
    #copies param file to data dir
    copy(Bio::Root::IO->catfile($exe_dir,"projector_parameters.txt"), $data_dir);

    my $call = 
        "$exe "
        ."1 "
        ."$prefix "
        ."30000000 "
        ."$fasta_file "
        ."$data_dir/ "
        ."$data_dir/ "
        ."projector_parameters.txt "
        ."$gtf_file "
        ."> "
        ."$log";
    system($call);
    print STDERR "Starting... $prefix\n";
} while(@exec_files);

1;

