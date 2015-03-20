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

say STDERR "## ::3:: MSA OF PROTEINS"                                                                                        ;
$alignment->runMuscle;

say STDERR "## ::4:: BUILD CONTIG MSA";
$alignment->buildContigMSA;

say STDERR "## ::4:: MAPPING CONTIGS TO REF PROT MSA"                                                                        ;

my $Finalout=Bio::SeqIO->new(-file=>">$out.msa",-format=>"fasta")                                                           ;
$Finalout->width(10000);
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
