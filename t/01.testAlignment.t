#!/usr/bin/env perl

use strict;
use warnings;

use FindBin qw/$Bin/;
use local::lib "$Bin/../local";
use Modern::Perl '2015';
use experimental qw/postderef signatures/;
use PASS::Alignment;
use Test::More;
use FindBin;
use Data::Dumper;
my $module = 'PASS::Alignment';

say "Testing PASS::Alignment package...";
use_ok($module)
    or BAIL_OUT("Does $module.pm compile?  Does it end with 1; ?");

say "Initialise object";
my $alignment = PASS::Alignment->new(
    refseqFasta     =>  "$FindBin::Bin/data/refSeqProtDB/ko\:K00001", #reference sequences
    alignmentFiles  =>  [glob("$FindBin::Bin/data/pAss01/K00001/alignment-*")],
    outputPrefix    =>  "$FindBin::Bin/data/pAss",
)                                                                              ;

say "Test methods:\n\t".'$obj->storeRefseq';
$alignment->storeRefseq();
ok scalar keys $alignment->{refseq}->%* == 772,
    "Loaded the correct number of contigs";

say "\t",'$obj->assignContig2ref';
$alignment->assignContig2ref();
my @missing;
for my $contigKey (keys $alignment->{contigs}->%*)
{
    for my $ref ($alignment->{contigs}->{$contigKey}->{parentREF})
    {
        if ($ref eq "")
        {
            push @missing, $contigKey;
        }
    }
}
is scalar @missing, 0, "Contigs matched to best aligned reference sequence";


say "\t",'$obj->runMuscle';
$alignment->runMuscle("data/pAss03/K00001.temp.ref.msa");
my @protMSAPositions = sort {$a <=> $b} values $alignment->{refseq}{'ref|NP_711797.1|'}{map}->%*;

    my $msa = <DATA>; chomp $msa;
    my ($lastLetter) = $msa =~ m/([A-Z])-+$/;
    #zero indexed
    my $position = $-[0];
    $position++;
    my @aaloc = sort {$a <=> $b} keys $alignment->{refseq}{'ref|NP_711797.1|'}{map}->%*;
    $msa =~ s/-//g;

is scalar @aaloc         , length $msa , "Load in protein MSA sequences";
is $protMSAPositions[-1] , $position   , "last letter is the same MSA loc";
my $lastLetterObj = substr $alignment->{refseq}{'ref|NP_711797.1|'}{seq}, -1;
is $lastLetterObj        , $lastLetter , "last letter is the same alphabet";


##buildContigMSA
#what if the reference sequence was not >ref|XX_XXXXXX|
#it seems like all the sequences at least for K00001 have the above structure

#Part1:
say "Refseq of interest:\t\t\tref|NP_711797.1|";
say "Contig: $_ Parent: \t\t$alignment->{contigs}{$_}{parentREF}" for qw/contig00025 contig01042 contig01169/;

is $alignment->{refseq}{'ref|NP_711797.1|'}{map}{107}, 246, "the protein MSA is stored properly";
$alignment->buildContigMSA();

is $alignment->{contigs}{'contig01042'}{globalCoordinates}{246}, 'GCG', 'contig01042 match';
is $alignment->{contigs}{'contig01169'}{globalCoordinates}{466}, 'GTT', 'contig01169 match';
is $alignment->{contigs}{'contig00025'}{globalCoordinates}{278}, 'ggt', 'contig00025 match';
#use 'ref|NP_711797.1|' to test again.
#data/pAss01/K00001/alignment-Leptospira_interrogans-ref_NP_711797.1__alcohol_dehydrogenase__Leptospira_interrogans_serovar_Lai_str.-01029.fasta
#data/pAss01/K00001/alignment-cellular_organisms-ref_NP_711797.1__alcohol_dehydrogenase__Leptospira_interrogans_serovar_Lai_str.-01553.fasta

my $Finalout=Bio::SeqIO->new(-file=>">K00001.contig.msa",-format=>"fasta");
#$Finalout->width(10000);
my @coi = qw/contig00025 contig01042 contig01169/;
#say Dumper $alignment->{contigs}{'contig00025'};
foreach my $contigID (sort keys $alignment->{contigs}->%*)
{
    if($contigID ~~ @coi){
    my @contigSeq;
    my $contigObj = Bio::Seq->new(
        -display_id => $contigID,
        -alphabet =>'dna',
    );
    say $contigID;
    my $parent      = $alignment->{contigs}{$contigID}{parentREF};
    my %contigHash  = $alignment->{contigs}{$contigID}{globalCoordinates}->%*;
#
    for my $aaXaa (sort {$a <=> $b } keys $alignment->{pos}->%*)
    {
#        say $aaXaa; #this is the positions for the whole protein msa
        my $codon = exists $contigHash{$aaXaa} ? $contigHash{$aaXaa} : '---';
        push @contigSeq,$codon;
    }
    my $fullsequence = join '', @contigSeq;
    #say $fullsequence;
    $contigObj->seq($fullsequence);
    $contigObj->desc($parent);
    $Finalout->write_seq($contigObj);
    }
}


##grepRefSeqID

##minmax

done_testing();

__DATA__
------------------------------------------------MTSLRIFKQVPRLLF----------------GFNTI-------------DRINELLPK-KNNGDYYIFIID--DVHQKGTIHSRLKHASED-----------------------------MIEWFPASVKEPSTLQ-ID----NLRD--------------KFMQARDHK-----L-----PKAIVGIGGG-STMDVAKA-LSVMMCNEGS---------ASQYQ----GWDLV-PNPGI------YKIGIP-------TVAGSGAEASRTAVLMGKERKFGINS-----------DHS-MFDAIILDS--SLIKNVPIA-QRFYS------GMDCYIHCVESLQ----------------------------GTMINELA-----KGNASKALELC------------------EKVFLSDGDDD---MLLT-ASYMGG---VSIV--N-----------------------SEVGV-CHALSYG------------LSLELGYRHG-------------------------------FANCVAFNVLDEY----Y--------------------------GP-----WVDR---------------FREML-------KIH--------KIELPK----NV---------CRSLDEAA--MQR------MVNMTLK----------MER------------------PLTNA--LGENWKDKMTP----NKIVSLYERM-------------------------------
>ref|NP_711797.1|
