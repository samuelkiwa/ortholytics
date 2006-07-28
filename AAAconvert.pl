use strict;
use DBI;
use IO::Handle;
use File::Find;
use Bio::Perl;
use Bio::SeqIO;
use Bio::Tools::GFF;

use constant SEARCH_DIRS => ".", "../seqs/AAA"; #Path to the files

my ($cnt, $cnt1, $line, $flag, $seq, $introns, @array, @array1);
my ($in_gen, $in_cds, $in_prot, $gen, $cds, $prot, @files);
my ($sp, $cod, @desc, @gen_lim, @gen_list);

#gets all the FASTA ans GFF files in the $path directories
find( sub{ m/\.genome\.fasta/ and unshift(@files,$File::Find::name);}, SEARCH_DIRS, );
$cnt1 = @files;
for($cnt=0; $cnt<$cnt1; $cnt++)
{
    @desc = split(m/\.|\//, $files[$cnt]);
    find( sub{ m/$desc[-3]\.gff$/ and push(@files,$File::Find::name);}, SEARCH_DIRS, );
    find( sub{ m/$desc[-3]\.orthologs\.cdna\.fasta$/ and push(@files,$File::Find::name);}, SEARCH_DIRS, );
    find( sub{ m/$desc[-3]\.orthologs\.translation\.fasta$/ and push(@files,$File::Find::name);}, SEARCH_DIRS, );
}
@array = sort {(split "/", $a)[-1] cmp (split "/", $b)[-1]} @files;
@files = @array;

print join("\n", @files);print("\n");

#Checks if exists all 4 files for every specie
$cnt = @files;
if($cnt1 == 0 || $cnt1*4 != $cnt){print("Error: files missing!!\n");exit(0);}
print (($cnt/4)." species found\n");

for($cnt=0; $files[$cnt] ne undef; $cnt+=4)
{
    #Opens the files
    $in_cds = new Bio::SeqIO(-file => $files[$cnt+2], -format => 'Fasta');
    $in_prot = new Bio::SeqIO(-file => $files[$cnt+3], -format => 'Fasta');
       
    #Gets species info
    @desc = split(m/\.|\//, $files[$cnt]);
    $cod = $sp = substr($desc[-3],0,3);
    $sp = "d".$sp;
    $cod =~ tr/a-z/A-Z/;

    #Opens the output files
    open(INGFF, "<$files[$cnt+1]");
    open(OUTGENE, ">$sp-all-gene-r0.0.fasta");
    open(OUTRNA, ">$sp-all-transcript-r0.0.fasta");
    open(OUTCDS, ">$sp-all-CDS-r0.0.fasta");
    open(OUTPROT, ">$sp-all-translation-r0.0.fasta");
    
    #Cycles through the CDS and Prot files
GENE: while(($cds = $in_cds->next_seq()) && ($prot = $in_prot->next_seq()))
    {
        #Gets gene info
        @desc = split(m/:/, $cds->desc);
        $flag=0;
        $introns = "";
        seek(*INGFF, 0, 0);
        $in_gen = new Bio::SeqIO(-file => $files[$cnt], -format => 'Fasta');
            
        #Gets gene exon boundaries from GFF
        while (!eof(INGFF))
        {
	          $line = readline(INGFF);
            @array = split(m/\t/, $line);
            @array1 = split(m/"/, $array[8]);
            
            if($array[2] eq "orthology_gene" && $array1[1] eq substr($cds->id,0,-3))
            {
                next GENE if($array[0] ne $desc[0]);
                $gen_lim[0] = $array[3];
                $gen_lim[1] = $array[4];
            }
            
            elsif($array[2] eq "orthology_exon" && $array1[1] eq $cds->id)
       	    {
                next GENE if($array[0] ne $desc[0]);
                if($desc[2] eq "+")
                {$introns .= $array[3]."..".$array[4].",";}
                else
                {$introns = $array[3]."..".$array[4].",".$introns;}
                $flag+=1;
            }
            elsif($flag!=0){chop($introns);last;}
        }

        #gets gene sequence from chromosome
        while($gen = $in_gen->next_seq())
        {
            if($gen->id eq $desc[0])
            {
                print($gen->id." len=".$gen->length." / ".$cod.substr($cds->id,2,-3)." ".$gen_lim[0]." - ".$gen_lim[1]."\n");

                if(($desc[2] eq "+") && ($gen_lim[1]+3 <= $gen->length) && ($seq = $gen->subseq($gen_lim[1],$gen_lim[1]+3) =~ m/TGA|TAA|TAG/gi))
                {
                    $seq = $gen->subseq($gen_lim[0],$gen_lim[1]+3);
                    
                    $gen_lim[-1] += 3;
                    $desc[1] = join("\.\.", @gen_lim);
                    last;
                }
                elsif(($desc[2] eq "-") && ($gen_lim[0] > 3) && ($seq = $gen->subseq($gen_lim[0]-3,$gen_lim[0]) =~ m/TCA|TTA|CTA/gi))
                {
                    $seq = $gen->subseq($gen_lim[0]-3,$gen_lim[1]);

                    $gen_lim[0] -= 3;
                    $desc[1] = join("\.\.", @gen_lim);
                    last;
                }
                else
                {
                    $seq = $gen->subseq($gen_lim[0],$gen_lim[1]);

                    $desc[1] = join("\.\.", @gen_lim);
                    last;
                }
            }
        }
        
        if($desc[2] eq "-"){$seq = reverse_complement_as_string($seq);}

        $line = $cod.substr($cds->id,2,-3);
        if(!grep(/^$line$/, @gen_list))
        {
            print(OUTGENE ">".$line." type=gene; loc=$desc[0]:".(($desc[2] eq "-")?("complement("):(""))."$desc[1]".(($desc[2] eq "-")?(")"):(""))."; id=".$line."; name=".$line."; release=r0.0; species=$sp; len=".length($seq)."\n");
            print(OUTGENE $seq."\n\n");
            push(@gen_list, $line);
        }
        print(OUTCDS ">".$cod.substr($cds->id,2,)." type=CDS; loc=$desc[0]:".(($desc[2] eq "-")?("complement("):(($flag!=0)?("join("):("")))."$introns".(($desc[2] eq "-" || $flag!=0)?(")"):(""))."; id=".$cod.substr($cds->id,2,)."; name=".$cod.substr($cds->id,2,)."; release=r0.0; species=$sp; len=".$cds->length."\n");
        print(OUTCDS $cds->seq."\n\n");
        print(OUTPROT ">".$cod.substr($prot->id,2,)." type=protein; loc=$desc[0]:".(($desc[2] eq "-")?("complement("):(($flag!=0)?("join("):("")))."$introns".(($desc[2] eq "-" || $flag!=0)?(")"):(""))."; id=".$cod.substr($prot->id,2,)."; name=".$cod.substr($prot->id,2,)."; release=r0.0; species=$sp; len=".$prot->length."\n");
        print(OUTPROT $prot->seq."\n\n");

        if(substr($seq,-3,3) =~ m/TGA|TAA|TAG/i)
        {
            @array = split(m/\.\.|,/, $introns);

            if($desc[2] eq "+"){$array[-1] += 3;}
            else{$array[0] -= 3;}
            
            $introns = "";
            for($cnt1=0; $array[$cnt1] ne undef; $cnt1+=2)
            {$introns .= $array[$cnt1]."..".$array[$cnt1+1].",";}
            chop($introns);
            
            $seq = $cds->seq . substr($seq,-3,3);
        }
        else{$seq = $cds->seq;}

        print(OUTRNA ">".$line."-R".substr($cds->id,-1,1)." type=transcript; loc=$desc[0]:".(($desc[2] eq "-")?("complement("):(($flag!=0)?("join("):("")))."$introns".(($desc[2] eq "-" || $flag!=0)?(")"):(""))."; id=".$line."-R".substr($cds->id,-1,1)."; name=".$line."-R".substr($cds->id,-1,1)."; release=r0.0; species=$sp; len=".length($seq)."\n");
        print(OUTRNA $seq."\n\n");
    }
        
    close OUTGENE;
    close OUTRNA;
    close OUTCDS;
    close OUTPROT;
    close INGFF;
}

print ("OK;");
exit(0);

