#!/usr/bin/env perl

use lib "./lib";
use v5.20                                                                                                                   ;
use feature 'signatures', 'postderef';
no warnings 'experimental';
use autodie                                                                                                                 ;
use Bio::SeqIO                                                                                                              ;
use MSA::Alignment                                                                                                          ;
use Data::Dumper;

use Statistics::Basic qw(mean stddev)                                                                                       ;
use Set::Intersection                                                                                                       ;

my ($dbfastaProt, $msadir, $out) = @ARGV                                                                                    ;
die "usage: $0 db.fasta msa.dir out.msa\n" unless $#ARGV == 2 && -f $dbfastaProt && -d $msadir                              ;

$msadir =~ s|//|/|g;$msadir =~ s|/$||; #deal with weird path names

my %refseq;        #stores reference sequences
my $temp = "$out".rand().time()                                                                                             ;

#Part:1###########################################
print STDERR "# ::1:: READ REFERENCE PROTEIN SEQUENCES"                                                                              ;
##################################################
    my $in = Bio::SeqIO->new(-file=>$dbfastaProt, -format=>"fasta")                                                         ;
    while(my $seqObj = $in->next_seq)                                                                                       {
        my $id = $seqObj->display_id;
        $id =~ m/(?<refID>ref\|\S+?\|)/; #MEGAN only keeps this
        $refseq{$+{refID}} = {seq =>  uc $seqObj->seq, length =>  $seqObj->length}                                 };
    #min and max referenceSequence lengths
    say STDERR "!!DONE";

    my @len = map {$refseq{$_}->{length}} keys %refseq                                                                           ;
    my $len = {min => mean(@len) - 2 * stddev(@len), max => mean(@len) + 2 * stddev(@len)                                   };
    say "#\tMIN ref length: $len->{min}\n#\tMAX ref length: $len->{max}\n";

#Part:2###########################################
say STDERR "# ::2:: ASSIGN BEST MATCH"                                                                                             ;
##################################################
    my $good_hit;   #assigns contigs to best refseq match
    #Criteria:
    #       1. length must be within 2 SDs of the mean length left and right (will no)
    assignContig2ref($_) for glob "$msadir/alignment-*"                                                                     ;
say STDERR "#\t!!DONE\n";

#Part:3###########################################
say STDERR "# ::3:: MSA OF PROTEINS"                                                                                               ;
print STDERR "#\tWRITING REFSEQS"                                                                                           ;
##################################################
    my $MSAout = Bio::SeqIO->new(-file => ">$out.temp.ref.faa",-fasta => "fasta")                                              ;
    my %good_ref                                                                                                            ;
    for my $contigID (keys $good_hit->%*)                                                                                   {
        my $protRef = $good_hit->{$contigID}{parentREF}                                                                     ;
        next if $good_ref{$protRef}                                                                                         ;
        my $outputSeq = Bio::Seq->new(
            -display_id => $protRef,
            -seq        => $refseq{$protRef}->{seq},
            -alphabet   => 'protein')                                                                                       ;
        $MSAout->write_seq($outputSeq)                                                                                         ;
        $good_ref{$protRef}++                                                                                               ;}
    #########################
say STDERR "\t!!DONE";
    print STDERR "#\tGenerate MSA with Muscle
##################################################"                                                                                               ;
    #########################

    system "muscle -in $out.temp.ref.faa -out $out.temp.ref.msa"                                                            ;
    #########################
    say STDERR "##################################################\n";
print STDERR "#\tMAPPING"                                                                                               ;
    #########################
    my $map                                                                                                                 ;
    my $in=Bio::SeqIO->new(-file=>"$out.temp.ref.msa", -format=>"fasta")                                                    ;
    while(my $seqObj = $in->next_seq)                                                                                       {
            my $seq = $seqObj->seq                                                                                          ;
            my $id = $seqObj->display_id                                                                                    ;
            my $aaIndex++                                                                                                   ;
            while($seq =~ m/[^-]/g)                                                                                         {
                $aaIndex++;   #for each aa
                my $aaLOC = $-[0] + 1                                                                                       ;
                $map->{$id}{$aaIndex} = $aaLOC                                                                                ;}}
say STDERR "\t!!DONE";

#Part:4###########################################
say STDERR "# ::4:: MAPPING CONTIGS TO REF PROT MSA"                                                                                 ;
##################################################
my $Finalout=Bio::SeqIO->new(-file=>">$out.msa",-format=>"fasta");
#my @reference = map loopMSA($_, $out) for glob "$msadir/alignment-*"                                                                            ;

#for my $id (keys %hit)                                                         {
    #my @sequence                                                               ;
    #for my $pos(sort {$a <=> $b} keys %pos)                                    {
        #$codon = $hit{$id}->{$pos}                                             ;
        #push @sequence, $codon ? $condon : '---'                               ;}
        #my $seqObj= Bio::Seq->new(-display_id =>$id, -seq=>join '', @sequence) ;}

#Because the same contigs can be aligned to different reference protein sequences
#So we only assign the contig to ref prot sequence if the the alignment requires the least amount of gaps to be introduced.
sub assignContig2ref ($alignment)                                                           {
    my $msaFile = (split /\//, $alignment)[-1]                                              ;
    say STDERR "#\tprocessing $msaFile"                                                        ;
    my $in=Bio::SeqIO->new(-file=>$alignment,-format=>'fasta')                              ;
## PROTEIN SEQUENCE
    my $refObj = $in->next_seq                                                              ;
    my $id = $refObj->display_id                                                            ;
    $id =~ m/(?<refID>ref\|\S+?\|)/; #MEGAN only keeps this
    #skip if the sequence is too short or too long ??WHY??

    return if $refseq{$+{refID}}->{length} < $len->{min} || $refseq{$+{refID}}->{length} > $len->{max}  ;
    my $protRef = $+{refID}                                                                       ;
## PARSE CONTIG MSA (DNA)
        while(my $seqObj = $in->next_seq)                                                   {
            my $contigID  =  $seqObj->display_id                                            ;
            my $sequence  =  $seqObj->seq                                                   ;
            my $numGaps   =  ($sequence =~ s/-//g)                                          ;
            my $len       =  $seqObj->length                                                ;
            $len    -=  $numGaps/1e6                                                        ;
            my $isBetterMatch = $len > $good_hit->{$contigID}{length}                       ;
            $good_hit->{$contigID} = {parentREF => $protRef, length => $len} if $isBetterMatch }}

sub loopMSA ($alignmentFile, $out)                                                                                          {
    #Function processes a single alignment file from MEGAN which is split into the reference seq and contigs
    my $in=Bio::SeqIO(-file   => $alignmentFile,-format => "fasta")                                                         ;
    #Step1: Reference Protein sequence
    my $seqObj = $in->next_seq                                                                                                 ;
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
        )                                                                                               ;
        $alignment->processRef                                                                                              ;
        #Step2: NT Contigs
        while($seqObj = $in->next_seq)                                                                                      {
            my $isComplement= $good_hit->{$seqObj->display_id}{id} eq $ref                                                  ;
            $isComplement ? $alignment->processContig($seqObj->seq, $seqObj->display_id) : next                                             ;}
        $alignment->dump;
    }else{
        return                                                                                                            }}
