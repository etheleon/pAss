use strict;
use warnings;

use Modern::Perl '2015';
use experimental qw/postderef signatures/;
use PASS::Alignment;
use Test::More;
use FindBin;
use Data::Dumper;
my $module = 'PASS::Alignment';

use_ok($module)
    or BAIL_OUT("Does $module.pm compile?  Does it end with 1; ?");

#Check functions are present
can_ok($module, 'storeRefseq')
    or BAIL_OUT("Missing package $module; or missing sub storeRefseq()");
can_ok($module, 'assignContig2ref')
    or BAIL_OUT("Missing package $module; or missing sub assignContig2ref()");
can_ok($module, 'runMuscle')
    or BAIL_OUT("Missing package $module; or missing sub runMuscle()");
can_ok($module, 'buildContigMSA')
    or BAIL_OUT("Missing package $module; or missing sub buildContigMSA()");
can_ok($module, 'grepRefSeqID')
    or BAIL_OUT("Missing package $module; or missing sub grepRefSeqID()");
can_ok($module, 'minmax')
    or BAIL_OUT("Missing package $module; or missing sub minmax()");


#Build an alignment OBJ
my $alignment = PASS::Alignment->new(
    refseqFasta     =>  "$FindBin::Bin/data/refSeqProtDB/ko\:K00001", #reference sequences
    alignmentFiles  =>  [glob("$FindBin::Bin/data/pAss01/K00001/alignment-*")],
    outputPrefix    =>  "$FindBin::Bin/data/pAss",
)                                                                              ;

#TEST METHODS
## storeRefseq
$alignment->storeRefseq();

ok scalar keys $alignment->{refseq}->%* == 772,
    "Loading: read allContigs In";

## assignContig2ref
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

is scalar @missing, 0, "Everyone is partnered";

$alignment->runMuscle("data/pAss03/K00001.temp.ref.msa");

my @protMSAPositions = sort {$a <=> $b} values $alignment->{refseq}{'ref|NP_711797.1|'}{map}->%*;


    my $msa = <DATA>; chomp $msa;

    my ($lastLetter) = $msa =~ m/([A-Z])-+$/;
#zero indexed
    my $position = $-[0];
    $position++;

    my @aaloc = sort {$a <=> $b} keys $alignment->{refseq}{'ref|NP_711797.1|'}{map}->%*;
    $msa =~ s/-//g;

is scalar @aaloc, length $msa, "protein sequences are the same length";
is $protMSAPositions[-1], $position, "last letter is the same MSA loc";
my $lastLetterObj = substr $alignment->{refseq}{'ref|NP_711797.1|'}{seq}, -1;
is $lastLetterObj, $lastLetter, "last letter is the same alphabet";


##buildContigMSA
#what if the reference sequence was not >ref|XX_XXXXXX|
#it seems like all the sequences at least for K00001 have the above structure

#Part1:
my $fullAlignment = "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------S--L--D--A--L--G--S--P--V--T--C--V--A--S--V--N--G--L--R--R--R--G--R--H--V--Q--V--G--L--L--P--S--S--T--G--S-----T--P--V--P--M--A--R--V--I--G--L--E--L--E--V--L--G--S--H--G--M--A--A--H--A--Y--P--P--M--L--E--L--";
my $fullRefseqSeq = "MRAVVFERFGEPAEVREVADPEPSEHGVVVRVEATGLCRSDWHGWVGHDPDITLPHVPGHELAGVVEAVGTRVSGWRPGDRVTVPFVCACGSCAACAAGDQQVCERQAQPGFTHWGSFAQYVALEHADVNLVAMPDELSFGTAAALGCRFATAFRAVVAQGRVAAGEWVAVHGCGGVGLSAVMIAAASGARVVAVDVSPQALDLARKFGAAASVDASRVEDTAAAVRELTGGGAHLSLDALGSPVTCVASVNGLRRRGRHVQVGLLPSSTGSTPVPMARVIGLELEVLGSHGMAAHAYPPMLELVRAGVLRPDLLVTSAIPLDTAPIALAAMGTAPGAGVTVIEPWR";
my $obj;

my $offset = 0;
while($fullAlignment =~ m/(?<continuousAA>(?<codon>[^-]--)+)/g)
{
    #will only match X-- patterns so wont work for ----
    my $start           =  $-[0];  #zero indexed position where the continuous stretch of amino acid position begin
    say "start: ", $start;
    say "start loc:", substr($fullAlignment, eval{$start - 3});
    my $spacedSeq       =  $+{continuousAA};
    (my $continuousSeq  =  $spacedSeq) =~ s/-//g;
    my $index           =  index(substr($fullRefseqSeq, $offset), $continuousSeq);  #incomplete alignment with the reference sequernce
    my $aa              =  0;
    while($spacedSeq =~ m/[^-]/g)
    {
    $aa++;                          #move stepwise from 1 aa to the next aa
    my $ntLoc = $start + $-[0] + 1; #plus one cause ltr we will do substr
    $obj->{ntMSALoc}{$ntLoc} = $aa + $offset + $index; #AA position
    }
    $offset += length $continuousSeq;
}
say "hello";
#for (sort {$a <=> $b} keys $obj->{ntMSALoc}->%*){say "$_\t$obj->{ntMSALoc}{$_}}


#use 'ref|NP_711797.1|' to test again.
data/pAss01/K00001/alignment-Leptospira_interrogans-ref_NP_711797.1__alcohol_dehydrogenase__Leptospira_interrogans_serovar_Lai_str.-01029.fasta
data/pAss01/K00001/alignment-cellular_organisms-ref_NP_711797.1__alcohol_dehydrogenase__Leptospira_interrogans_serovar_Lai_str.-01553.fasta


##grepRefSeqID

##minmax

done_testing();

__DATA__
------------------------------------------------MTSLRIFKQVPRLLF----------------GFNTI-------------DRINELLPK-KNNGDYYIFIID--DVHQKGTIHSRLKHASED-----------------------------MIEWFPASVKEPSTLQ-ID----NLRD--------------KFMQARDHK-----L-----PKAIVGIGGG-STMDVAKA-LSVMMCNEGS---------ASQYQ----GWDLV-PNPGI------YKIGIP-------TVAGSGAEASRTAVLMGKERKFGINS-----------DHS-MFDAIILDS--SLIKNVPIA-QRFYS------GMDCYIHCVESLQ----------------------------GTMINELA-----KGNASKALELC------------------EKVFLSDGDDD---MLLT-ASYMGG---VSIV--N-----------------------SEVGV-CHALSYG------------LSLELGYRHG-------------------------------FANCVAFNVLDEY----Y--------------------------GP-----WVDR---------------FREML-------KIH--------KIELPK----NV---------CRSLDEAA--MQR------MVNMTLK----------MER------------------PLTNA--LGENWKDKMTP----NKIVSLYERM-------------------------------
>ref|NP_711797.1|
