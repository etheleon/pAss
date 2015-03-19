#!/usr/bin/env perl
use Modern::Perl '2015'                                                                                                     ;
use experimental qw/signatures postderef/                                                                                       ;
use autodie                                                                                                                 ;
use PASS::Alignment                                                                                                          ;
use Data::Dumper                                                                                                            ;

my ($refseqFasta, $msadir, $out) = @ARGV                                                                                    ;
die "usage: $0 db.fasta msa.dir out.msa\n" unless $#ARGV == 2 && -f $refseqFasta && -d $msadir                              ;
$msadir =~ s|//|/|g                                                                                                         ;
$msadir =~ s|/$||                                                                                                           ;

my $alignment = PASS::Alignment->new(
    refseqFasta     =>  $refseqFasta, #reference sequences
    alignmentFiles  =>  [glob("$msadir/alignment-*")],
    outputPrefix    =>  $out,
);

say STDERR "## ::1:: STORING REFERENCE PROTEIN SEQUENCES"                                                                     ;
$alignment->storeRefseq;

say STDERR "## ::2:: ASSIGNING CONTIGS TO BEST REFSEQ ALIGNMENT"                                                                                      ;
$alignment->assignContig2ref;
#say Dumper $alignment->{contigs}{'contig00049'};

say STDERR "## ::3:: MSA OF PROTEINS"                                                                                        ;
$alignment->runMuscle;

say STDERR "## ::4:: BUILD CONTIG MSA";
$alignment->buildContigMSA;

say STDERR "## ::4:: MAPPING CONTIGS TO REF PROT MSA"                                                                        ;

my $Finalout=Bio::SeqIO->new(-file=>">$out.msa",-format=>"fasta")                                                           ;

my %position = $alignment->{pos}->%*;
my $msalength = scalar keys %position;
say $msalength;
exit 1;

$Finalout->width($msalength);
#say Dumper $alignment->{pos};
foreach my $contigID (keys $alignment->{contigs}->%*)
{
    my @contigSeq;
    my $contigObj = Bio::Seq->new(
        -display_id => $contigID,
        -alphabet =>'dna',
    );
    my $parent = $alignment->{contigs}{$contigID}{parentREF};
    my %contigHash = $alignment->{contigs}{$contigID}{globalCoordinates}->%*;

    for my $aaXaa(sort {$a <=> $b } keys $alignment->{pos}->%*)
    {
        my $codon = exists $contigHash{$aaXaa} ? $contigHash{$aaXaa} : '---';
        push @contigSeq,$codon;
    }
    my $fullsequence = join '', @contigSeq;
    $contigObj->seq($fullsequence);
    $contigObj->desc($parent);
    $Finalout->write_seq($contigObj);
}
#my @reference = map loopMSA($_, $out) for glob "$msadir/alignment-*"                                                       ;

#for my $id (keys %hit)                                                                                                     {
#my @sequence                                                                                                           ;
#for my $pos(sort {$a <=> $b} keys %pos)                                                                                {
#$codon = $hit{$id}->{$pos}                                                                                         ;
#push @sequence, $codon ? $condon : '---'                                                                           ;}
#my $seqObj= Bio::Seq->new(-display_id =>$id, -seq=>join '', @sequence)                                             ;}

#Because the same contigs can be aligned to different reference protein sequences
#So we only assign the contig to ref prot sequence if the the alignment requires the least amount of gaps to be introduced.

#sub loopMSA ($alignmentFile, $out)                                                                                          {
#    #Function processes a single alignment file from MEGAN which is split into the reference seq and contigs
#    my $in=Bio::SeqIO(-file   => $alignmentFile,-format => "fasta")                                                         ;
#    #Step1: Reference Protein sequence
#    my $seqObj = $in->next_seq                                                                                              ;
#    my $ref             =  $seqObj->display_id                                                                              ;
#    my $isReferenceSeq  =  ($ref =~ m/^ref\|/)                                                                              ;
#    my $isGood          =  exists $good_ref{$ref}                                                                           ;
#    if($isReferenceSeq && $isGood)                                                                                          {
#        #BUILD globalCoordinates
#        my $alignment = MSA::Alignment->new(
#            sequence => $seqObj->seq,
#            id       => $ref,
#            output   =>$out,
#            map => $map,
#            #hit => $hit,
#        )                                                                                                                   ;
#        $alignment->processRef                                                                                              ;
#        #Step2: NT Contigs
#        while($seqObj = $in->next_seq)                                                                                      {
#            my $isComplement= $good_hit{$seqObj->display_id}->{id} eq $ref                                                  ;
#            $isComplement ? $alignment->processContig($seqObj->seq, $seqObj->display_id) : next                             ;}
#        $alignment->dump                                                                                                    ;
#    }else                                                                                                                   {
#        return                                                                                                              }}
