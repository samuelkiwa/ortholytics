# create_ribosomal_sets.PLS
#
# Cared for by Albert Vilella <>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

create_ribosomal_sets.PLS - Obtains CDS nucleotide sequence for
ribosomal entries from a rnt file (similar to ptt but constrained to
ribo-related products) and the corresponding genbank file. The CDS are
exported to a single file in fasta format.

=head1 SYNOPSIS

perl
create_ribosomal_sets.PLS
-rnt /path/to/file.rnt
-gbk /path/to/file.gbk
-o /path/to/outdir

=head1 DESCRIPTION

Describe the object here

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

my ($rntinput,$gbkinput,$outdir) = "";
my $out;

GetOptions(
	   'rntinput|rnt|r:s' => \$rntinput,
	   'gbkinput|gbk|g:s' => \$gbkinput,
	   'o|outdir:s'=> \$outdir,
          );

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $logfile_time = 
    sprintf ("%04d%02d%02d_%02d%02d%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec);

print "Reading informant rnt/ptt file...\n";
open(RNTFILE, $rntinput) or die "cannot open rnt/ptt file $rntinput: $!";
my $reset = 1;
my $count = 0;
my ($strand, $start, $end, $longname);
my %entries;
while (<RNTFILE>) {
    #           start    end    strand  len     gi       name
    if ($_ =~ /(\d+)\.\.(\d+)\s+([+-])\s+\d+\s+(\d+).+\t(.+)$/) {
        $start = sprintf("%08d", $1);
        $end   = sprintf("%08d", $2);
        $strand = $3;
        $longname = join ("_", split (/\s+/, $5));
        my $entry_name = "_" . $start . "_" . $end . "_" . $strand;
        $entries{$entry_name}{start} = $start;
        $entries{$entry_name}{end} = $end;
        $entries{$entry_name}{strand} = $strand;
        $entries{$entry_name}{longname} = $longname;
        $count++;
    }
}

close(RNTFILE);
print "Loaded $count entries...\n"; 
$count = 0;

my $gbkfile = new Bio::SeqIO(-file => $gbkinput,
                             -format => 'genbank');

unless (-d "$outdir") {
    mkdir "$outdir" or die "could not create directory $outdir: $!";
}

$out = Bio::SeqIO->new(-file=>Bio::Root::IO->catfile(">","$outdir","ribosomal.$logfile_time.fasta"),
                       -format => 'fasta');

while( my $seq = $gbkfile->next_seq ) {
    foreach my $entry (sort keys %entries) {
        my $entryseq;
        if ($entries{$entry}{strand} == "-") {
            $entryseq = $seq->revcom->subseq($entries{$entry}{start},$entries{$entry}{end});
        } else {
            $entryseq = $seq->subseq($entries{$entry}{start},$entries{$entry}{end});
        }
        my $length = length $entryseq;
        my $id = $entries{$entry}{longname};
        my $riboseq = new Bio::LocatableSeq(-seq => $entryseq,
                                            -id  => $id,
                                            -start => 1,
                                            -end   => $length
                                           );
        $out->write_seq($riboseq);
        my $pseq = $riboseq->translate(undef, undef, undef, 11, 'fullCDS', undef, undef);
        print "$pseq->seq\n\n";
    }
}

$out->close();

print "Done\n";

1;


