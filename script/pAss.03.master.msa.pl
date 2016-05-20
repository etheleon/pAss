#!/usr/bin/env perl

use FindBin qw/$Bin/;
use local::lib "$Bin/local/";
use Modern::Perl '2015'                                                        ;
use Pod::Usage;
use experimental qw/signatures postderef/                                      ;
use autodie                                                                    ;
use PASS::Alignment                                                            ;
#use Data::Dumper                                                              ;
use Getopt::Lucid qw/:all/                                                     ;

my @specs = (
    Param("ko|k"),
    Param("refseqFasta|f"),
    Param("msadir|m"),
    Param("out|o"),
    Switch("help|h")
)                                                                              ;

my $opt = Getopt::Lucid->getopt( \@specs )                                     ;
pod2usage(-verbose=>2) if $opt->get_help                                       ;
$opt->validate({'requires' => ['refseqFasta', 'msadir', 'out']})               ;

my $refseqFasta = $opt->get_refseqFasta                                        ;
my $msadir      = $opt->get_msadir                                             ;
my $out         = $opt->get_out                                                ;
system "mkdir -p $out" unless -d "$out";
my $ko = $opt->get_ko;

$out .= "/$ko";

my $alignment = PASS::Alignment->new(
    refseqFasta     =>  $refseqFasta, #reference sequences
    alignmentFiles  =>  [glob("$msadir/alignment-*")],
    outputPrefix    =>  $out,
)                                                                              ;

say STDERR "## ::1:: STORING REFERENCE PROTEIN SEQUENCES"                      ;
$alignment->storeRefseq                                                        ;

say STDERR "## ::2:: ASSIGNING CONTIGS TO BEST REFSEQ ALIGNMENT"               ;
$alignment->assignContig2ref                                                   ;

say STDERR "## ::3:: MSA OF PROTEINS"                                          ;
$alignment->runMuscle                                                          ;

say STDERR "## ::4:: BUILD CONTIG MSA"                                         ;
$alignment->buildContigMSA                                                     ;

say STDERR "## ::4:: MAPPING CONTIGS TO REF PROT MSA"                          ;

my $Finalout=Bio::SeqIO->new(-file=>">$out.msa",-format=>"fasta")              ;
$Finalout->width(10000)                                                        ;
foreach my $contigID (sort keys $alignment->{contigs}->%*)
{
    my @contigSeq                                                              ;
    my $contigObj = Bio::Seq->new(
        -display_id => $contigID,
        -alphabet =>'dna',
    )                                                                          ;
    my $parent = $alignment->{contigs}{$contigID}{parentREF}                   ;
    $parent .= " ".$alignment->{contigs}{$contigID}{direction}
    my %contigHash = $alignment->{contigs}{$contigID}{globalCoordinates}->%*   ;

    for my $aaXaa(sort {$a <=> $b } keys $alignment->{pos}->%*)
    {
        my $codon = exists $contigHash{$aaXaa} ? $contigHash{$aaXaa} : '---'   ;
        push @contigSeq,$codon                                                 ;
    }
    my $fullsequence = join '', @contigSeq                                     ;
    $contigObj->seq($fullsequence)                                             ;
    $contigObj->desc($parent)                                                  ;
    $Finalout->write_seq($contigObj)                                           ;
}

=pod

=head1 NAME
    maxDivMain -

=head1 OPTIONS

=over 4

=item --ko -k

    the ko

=item --refseqFasta -f

    Path to FASTA file storing refseq protein sequences from belonging to KO
    eg. /export2/home/uesu/db/konr/ko:K00001

=item --msadir -m

    Path to MEGAN's alignment output ie. from pAss.01. Kept within each folder
    eg. /export2/home/uesu/reDiamond/out/assm.0400/K00001

=item --out -o

    Output Directory:
    Inside the output directory there will be

=back

=cut

