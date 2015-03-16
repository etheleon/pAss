#!/usr/bin/env perl
use Modern::Perl '2015'                                                                                                     ;
use feature 'signatures', 'postderef'                                                                                       ;
use autodie                                                                                                                 ;
use Bio::SeqIO;
use PASS::Alignment                                                                                                          ;
use Statistics::Basic qw(mean stddev)                                                                                       ;
use Set::Intersection                                                                                                       ;
use Data::Dumper                                                                                                            ;

my ($dbfastaProt, $msadir, $out) = @ARGV                                                                                    ;
die "usage: $0 db.fasta msa.dir out.msa\n" unless $#ARGV == 2 && -f $dbfastaProt && -d $msadir                              ;

$msadir =~ s|//|/|g                                                                                                         ;
$msadir =~ s|/$||                                                                                                           ;
my $refseqRegex= qr/(?<refseqID>ref\|\S+?\|)/                                                                                     ;

say STDERR "## ::1:: STORE REFERENCE PROTEIN SEQUENCES"                                                                     ;
    my %refseq                                                                                                                  ;
    my $dbin = Bio::SeqIO->new(-file=>$dbfastaProt, -format=>"fasta")                                                       ;
    while(my $seqObj = $dbin->next_seq)                                                                                     {
        my $id = $seqObj->display_id                                                                                        ;
        $id =~ m/$refseqRegex/                                                                                     ;
        $refseq{$+{refseqID}} = {seq => $seqObj->seq, length => $seqObj->length}                                          };

say STDERR "## ::2:: ASSIGN CONTIGS TO BEST REFSEQ ALIGNMENT"                                                                                      ;
    my %good_hit;
    my $len = minmax();
    assignContig2ref($_) for glob "$msadir/alignment-*"                                                                     ;

say STDERR "## ::3:: MSA OF PROTEINS"                                                                                        ;
    say STDERR "#\t ::a:: WRITING REFSEQS"                                                                                           ;
        my $MSAout = Bio::SeqIO->new(-file => ">$out.temp.ref.faa", -fasta => "fasta")                                           ;
        my %good_ref                                                                                                            ;
        foreach my $contig (keys %good_hit)                                                                                   {
            my $protRef = $good_hit{$contig}->{parentREF}                                                                     ;
            next if exists $good_ref{$protRef}                                                                                         ;
            my $outputSeq = Bio::Seq->new(
                -display_id => $protRef,
                -seq        => $refseq{$protRef}->{seq},
                -alphabet   => 'protein')                                                                                       ;
            $MSAout->write_seq($outputSeq)                                                                                      ;
            $good_ref{$protRef}++                                                                                               ;}

    say STDERR "##\t ::b:: Generate MSA with Muscle";
    runMuscle();

    say STDERR "##\t ::c:: MAPPING"                                                                                               ;
        my $map                                                                                                                 ;
        my $in=Bio::SeqIO->new(-file=>"$out.temp.ref.msa", -format=>"fasta")                                                    ;
        while(my $seqObj = $in->next_seq)                                                                                       {
                my $seq = $seqObj->seq                                                                                          ;
                my $id = $seqObj->display_id                                                                                    ;
                my $aaIndex++                                                                                                   ;
                while($seq =~ m/[^-]/g)                                                                                         {
                    $aaIndex++;   #for each aa
                    my $aaLOC = $-[0] + 1                                                                                       ;
                    $map->{$id}{$aaIndex} = $aaLOC                                                                              ;}}

say STDERR "## ::4:: MAPPING CONTIGS TO REF PROT MSA"                                                                        ;
    my $Finalout=Bio::SeqIO->new(-file=>">$out.msa",-format=>"fasta")                                                           ;
#my @reference = map loopMSA($_, $out) for glob "$msadir/alignment-*"                                                       ;

#for my $id (keys %hit)                                                                                                     {
    #my @sequence                                                                                                           ;
    #for my $pos(sort {$a <=> $b} keys %pos)                                                                                {
        #$codon = $hit{$id}->{$pos}                                                                                         ;
        #push @sequence, $codon ? $condon : '---'                                                                           ;}
        #my $seqObj= Bio::Seq->new(-display_id =>$id, -seq=>join '', @sequence)                                             ;}

#Because the same contigs can be aligned to different reference protein sequences
#So we only assign the contig to ref prot sequence if the the alignment requires the least amount of gaps to be introduced.
sub assignContig2ref ($alignment)                                                                                           {
    my $msaFile = (split /\//, $alignment)[-1]                                                                              ;
    my $in=Bio::SeqIO->new(-file=>$alignment,-format=>'fasta')                                                              ;
## PROTEIN SEQUENCE
    my $refObj = $in->next_seq                                                                                              ;
    my $id = $refObj->display_id                                                                                            ;
    $id =~ m/$refseqRegex/; #MEGAN only keeps this
    say $refseq{$+{refseqID}}->{length};
    say join "\t", $len->{min}, $len->{max};
    my $isCorrectLength = $refseq{$+{refseqID}}->{length} < $len->{min} || $refseq{$+{refseqID}}->{length} > $len->{max};
    return if $isCorrectLength;
    my $protRef = $+{refseqID}                                                                                                 ;
## PARSE CONTIG MSA (DNA)
    while(my $seqObj = $in->next_seq)                                                                                   {
        my $contigID  =  $seqObj->display_id                                                                            ;
        my $sequence  =  $seqObj->seq                                                                                   ;
        my $numGaps   =  ($sequence =~ s/-//g)                                                                          ;
        my $len       =  $seqObj->length                                                                                ;
        $len    -=  $numGaps/1e6                                                                                        ;
        if(exists $good_hit{$contigID}){
            my $isBetterMatch = $len > $good_hit{$contigID}->{length}                                                       ;
            $good_hit{$contigID} = {parentREF => $protRef, length => $len} if $isBetterMatch
        }else{
            $good_hit{$contigID} = {parentREF => $protRef, length => $len}
        }
    }}

sub loopMSA ($alignmentFile, $out)                                                                                          {
    #Function processes a single alignment file from MEGAN which is split into the reference seq and contigs
    my $in=Bio::SeqIO(-file   => $alignmentFile,-format => "fasta")                                                         ;
    #Step1: Reference Protein sequence
    my $seqObj = $in->next_seq                                                                                              ;
    my $ref             =  $seqObj->display_id                                                                              ;
    my $isReferenceSeq  =  ($ref =~ m/^ref\|/)                                                                              ;
    my $isGood          =  exists $good_ref{$ref}                                                                           ;
    if($isReferenceSeq && $isGood)                                                                                          {
        #BUILD globalCoordinates
        my $alignment = MSA::Alignment->new(
            sequence => $seqObj->seq,
            id       => $ref,
            output   =>$out,
            map => $map,
            #hit => $hit,
        )                                                                                                                   ;
        $alignment->processRef                                                                                              ;
        #Step2: NT Contigs
        while($seqObj = $in->next_seq)                                                                                      {
            my $isComplement= $good_hit{$seqObj->display_id}->{id} eq $ref                                                  ;
            $isComplement ? $alignment->processContig($seqObj->seq, $seqObj->display_id) : next                             ;}
        $alignment->dump                                                                                                    ;
    }else                                                                                                                   {
        return                                                                                                              }}

sub minmax {
    #Criteria for selecting refseq sequences
    #       1. length must be within 2 SDs of the mean length left and right (will no)
    my @len = map {$refseq{$_}->{length}} keys %refseq                                                                      ;
    my $len = {min => mean(@len) - 2 * stddev(@len), max => mean(@len) + 2 * stddev(@len)                                   };
    say "#\tMIN ref length: $len->{min}\n#\tMAX ref length: $len->{max}\n"                                                  ;
    return $len;
    }

sub runMuscle{
print STDERR '#' x 50;
system "muscle -in $out.temp.ref.faa -out $out.temp.ref.msa"                                                            ;
print STDERR '#' x 50;
}

