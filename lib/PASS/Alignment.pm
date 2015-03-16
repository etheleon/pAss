#!/usr/bin/env perl

package MSA::Alignment;
use Moo;
use v5.20;
use feature 'signatures', 'postderef';
no warnings 'experimental';
use strict;
use namespace::clean;

sub processRef($self){
    #gives a reference keypair of locations key:position value:aa
    my $fullseq = $self->sequence;
    my $refseq= $self->refseq;
    my $last_end = 0;
    while($fullseq=~ m/(([^-]--)+)/g){  #will only match A--
        my $start = $-[0];  #LOC before the 1st amino acid
        my $fullseq = $1;   #sequence after the gaps
        (my $strippedSeq = $fullseq) =~ s/-//g;
        my $index = index(substr($refseq, $last_end), $strippedSeq);  #sometimes it isnt a full alignment with the reference sequernce, so i need to know where it begins

        my $aa = 0;
        while($fullseq =~ m/[^-]/g)
        { #move stepwise from 1 aa to the next aa
            $aa++;
            my $ntLoc = $start + $-[0] + 1; #plus one cause ltr we will do substr
            $self->{msa2ref}{$ntLoc} = $aa + $last_end + $index;
        }
        $last_end += length $strippedSeq;}};

sub processContig($self, $contigSequence, $id)                                                                                                {
    while($contigSequence =~ m/([^-]{3})/g)                                                                                                   {
        my $codon                                 =  $1                                                                                       ;
        my $msaLoc                                =  $-[0] + 1                                                                                ;
        my $aaLOC                                 =  $self->{msa2ref}->{$msaLoc}                                                                        ;
        my $this                                  =  $self->{map}->{$id}{$aaLOC}                                                                     ;
        $this ? $self->position($this) : $self->upPos(1/1e6)                                                                                  ;
        $self->pos->{$self->position}++                                                       ;    #this logs the number of sequences in that particular positon
        $self->hit->{$id}{$self->position}                         =  $codon                                                                             ;}};

sub upPos($self, $amount){
    my $i = $self->position;
    $i = $i + $amount;
    $self->position($i);
};

sub dump($self){
    return {
        pos=>$self->pos,
        hit=>$self->hit}
}

has msa2ref => (
    is => 'rw',
    handles => {},
);

has map => (
    is => 'rw',
    handles => {},
);

has hit => (
    is => 'rw',
    handles => {},
);

has pos =>(
    is => 'rw',
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

#use strict;
#use v5.20;
#use experimental 'postderef';
#use Data::Dumper;
#use Bio::Seq;
#use Bio::SeqIO;
#use MSA::Alignment;

#say "#Test";
#my $ref= Bio::Seq->new(
    #-seq=>"PALGLGTWKSSPQVVGQAVEQALDLGYRHLDCAAIYGNEAEIGATLANAFTKGVVKREELWITSKLWSNAHHPDAVLPALEKTLQDLGLDYLDLYLIHWPVVIQPDVGFPESGDQLLPFTPASLEGTWQALEKAVDLGLCHHIGVSNFSLKKLEMVLSMARIPPAVNQVELHPYLQQSDLLTFANSQNILLTAYSPLGSGDRPAAFQQAAEPKLLTDPVINGIAAEQGCSAAQVLLAWAIQRGTVTIPKSVNPERLEQNLRAADITLTDSEMAKIALLDRHYRYVSGDFWTMPGSPYTLQNLWDE",
    #-display_id=>"ref",
    #-alphabet=>'protein', -format=>"fasta");
#my $in = Bio::SeqIO->new(-file=>"example/msaFile", -format=>"fasta");
#my $head= $in->next_seq;
#my $alignmentObj = MSA::Alignment->new(
    #sequence=>$head->seq,
    #id=>$head->display_id,
    #refseq=>$ref->seq,
#);

#$alignmentObj->processRef;
#while(my $seqObj = $in->next_seq){
#$alignmentObj->processContig($seqObj->seq,$seqObj->display_id);
#}

#say Dumper($alignmentObj);
