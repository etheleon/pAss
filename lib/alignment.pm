#!/usr/bin/env perl

package MSA::Alignment;
use Moo;
use v5.20;
use experimental 'signatures';
use experimental 'postderef';
use strict;
use namespace::clean;

sub processRef($self){
    #gives a reference keypair of locations key:position value:aa
    my $fullseq = $self->sequence;
    my $refseq= $self->refseq;
    my $last_end = 0;
    my %msa2ref;
    while($fullseq=~ m/(([^-]--)+)/g){  #will only match A--
        my $start = $-[0];  #LOC before the 1st amino acid
        my $fullseq = $1;   #sequence after the gaps
        (my $strippedSeq = $fullseq) =~ s/-//g;
        my $index = index(substr($refseq, $last_end), $strippedSeq);  #sometimes it isnt a full alignment with the reference sequernce, so i need to know where it begins

        my $aa = 0;
        while($fullseq =~ m/[^-]/g){ #move stepwise from 1 aa to the next aa
            $aa++;
            my $ntLoc = $start + $-[0] + 1; #plus one cause ltr we will do substr
            $msa2ref{$ntLoc} = $aa + $last_end + $index;}
        $last_end += length $strippedSeq;}
    $self->msa2ref(\%msa2ref);
};

#sub processContig($self, $contigSequence, $id)                             {
    #while($contigSequence =~ m/([^-]{3})/g)                                {
        #my $codon                                 =  $1                    ;
        #my $msaLoc                                =  $-[0] + 1             ;
        #my $aaLOC                                 =  $msa2ref{$msaLoc}     ;
        #my $this                                  =  $map->{$ref}{$aaLOC}     ;

        #$this ? $self->position($this) : $self->upPos(1/1e6);
        #$pos{$pos}++                                                       ;
        #$hit{$id}->{$pos}                         =  $codon                ;}}



sub upPos($self, $amount){
    my $i = $self->position;
    $i+=$amount;
    $self->position($i);
};

has msa2ref => (
    is => 'rw',
    isa => sub {die "Not valide" unless scalar keys $_[0]->%* },
);

has position => (
    is  => 'rw',
    default=> 0,
);

## Required attributes
has refseq => (
    is => 'ro',
    isa => sub {die "No spaces allowed" unless $_[0] =~ m/[^\s]/},
    required=>1,
);

has sequence => (
    is => 'ro',
    isa => sub {die "No spaces allowed" unless $_[0] =~ m/[^\s]/},
    required=>1,
);

has id => (
    is  => 'ro',
    required=>1,
);
1;

use strict;
use v5.20;
use experimental 'postderef';
use Data::Dumper;
use Bio::Seq;
use Bio::SeqIO;
use MSA::Alignment;

say "#Test";
my $ref= Bio::Seq->new(
    -seq=>"PALGLGTWKSSPQVVGQAVEQALDLGYRHLDCAAIYGNEAEIGATLANAFTKGVVKREELWITSKLWSNAHHPDAVLPALEKTLQDLGLDYLDLYLIHWPVVIQPDVGFPESGDQLLPFTPASLEGTWQALEKAVDLGLCHHIGVSNFSLKKLEMVLSMARIPPAVNQVELHPYLQQSDLLTFANSQNILLTAYSPLGSGDRPAAFQQAAEPKLLTDPVINGIAAEQGCSAAQVLLAWAIQRGTVTIPKSVNPERLEQNLRAADITLTDSEMAKIALLDRHYRYVSGDFWTMPGSPYTLQNLWDE",
    -display_id=>"ref",
    -alphabet=>'protein', -format=>"fasta");
my $in = Bio::SeqIO->new(-file=>"example/msaFile", -format=>"fasta");
my $seqObj = $in->next_seq;
my $alignmentObj = MSA::Alignment->new(
    sequence=>$seqObj->seq,
    id=>$seqObj->display_id,
    refseq=>$ref->seq,
);
my $good_hit;
#say Dumper($alignmentObj);
$alignmentObj->processRef;
my %hash = $alignmentObj->msa2ref->%*;
say "$_\t$hash{$_}" for keys %hash;
