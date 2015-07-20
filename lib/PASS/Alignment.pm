#!/usr/bin/env perl

package PASS::Alignment;
use v5.20;
use Moo;
use experimental qw/signatures postderef/;
use namespace::clean;
use Bio::SeqIO;
use Statistics::Basic qw/mean stddev/;
use List::MoreUtils qw/uniq/;
use Data::Dumper;

#main methods
#1
sub storeRefseq($self){
    my $dbin = Bio::SeqIO->new(
        -file=>$self->refseqFasta,
        -format=>"fasta");
    while(my $seqObj = $dbin->next_seq)
    {
        my $refseqID = $self->grepRefSeqID($seqObj->display_id);
        $self->{refseq}{$refseqID} = {seq => $seqObj->seq, length => $seqObj->length}
    };
}
#2
sub assignContig2ref($self){
    $self->minmax;

    for my $alignmentFile ($self->alignmentFiles->@*)
    {
        my $msaFile = (split /\//, $alignmentFile)[-1];
        my $in = Bio::SeqIO->new(-file=>$alignmentFile,-format=>'fasta');

#       PROTEIN SEQUENCE
        my $refObj = $in->next_seq;
        my $refseqID = $self->grepRefSeqID($refObj->display_id);
        #Some KOs cannot proceed because
        #1. left out cause their reference sequences are all the same length; SD is low and the alignment length is not the length
        #2. Too little reference sequences
        #etc.
        my $isCorrectLength = $self->{refseq}{$refseqID}{length} > $self->min &&  $self->{refseq}{$refseqID}{length} < $self->max;
        #refseq sequences shd be trimmed because ltr we will pivot the choice of alignment based on whther the sequence is too long;
        $self->{refseq}{$refseqID}{safe} = $isCorrectLength ? 1 : 0;
        next unless $isCorrectLength;


# PARSE CONTIG MSA (DNA)
        while(my $seqObj = $in->next_seq){
            my $contigID  =  $seqObj->display_id;
            my $sequence = $seqObj->seq;
            my $gap = ($sequence =~ s/-//g);
            my $length = length $sequence;
            $length -= $gap / 1e6;
            #say "$alignmentFile\n$refseqID\t$length" if $contigID eq 'contig00049';
            my $contigDetails = {
                seq  =>  $seqObj->seq,
                len =>   $length,
                parentREF => $refseqID,
            };
            if(!exists $self->{contigs}{$contigID}){
                $self->{contigs}{$contigID} = $contigDetails;
            }else{
                my $isBetterMatch = $contigDetails->{len} > $self->{contigs}{$contigID}{len};
                if ($isBetterMatch){
                    #say "Found better match";
                    $self->{contigs}{$contigID} = $contigDetails
                }}}}}
#3
sub runMuscle($self){
    #muscle is heuristic, each run will give you a different output
    my $out = $self->outputPrefix;
    my $MSAout = Bio::SeqIO->new(-file => ">$out.temp.ref.faa", -fasta => "fasta");
    $MSAout->width(1000);
    my @protRef = uniq map {$self->{contigs}{$_}{parentREF}} keys $self->contigs->%*;
    say "number of reference sequences included: ", scalar @protRef;
    foreach (@protRef){
        my $outputSeq = Bio::Seq->new(
            -display_id => $_,
            -seq        => $self->{refseq}{$_}{seq},
            -alphabet   => 'protein');
        $MSAout->write_seq($outputSeq);
    }

    say     STDERR "##\t ::b:: Generate MSA with Muscle";
    print   STDERR '#' x 50;
    system  "muscle -in $out.temp.ref.faa -out $out.temp.ref.msa";
    print   STDERR '#' x 50;
    say     STDERR "";

    my $MUSCLEinput=Bio::SeqIO->new(-file=>"$out.temp.ref.msa", -format=>"fasta")                                                    ;
    while(my $seqObj = $MUSCLEinput->next_seq)                                                                                       {
        my $seq = $seqObj->seq                                                                                          ;
        my $refseqID = $seqObj->display_id                                                                                    ;
        my $aaIndex = 1                                                                                                   ;
        while($seq =~ m/[^-]/g)
        {
            #$-[0] is zero index
            my $msaLOC = $-[0] + 1                                                                                       ;
            $self->{refseq}{$refseqID}{map}{$aaIndex} = $msaLOC                                                                              ;
            $aaIndex++;
        }
    }
}
#4
sub buildContigMSA($self){
    #gives a reference keypair of locations key:position value:aa
    foreach my $alignmentFile ($self->alignmentFiles->@*){
        my $in=Bio::SeqIO->new(-file => $alignmentFile, -format => "fasta");
        my $seqObj = $in->next_seq;
        my $refseqID             =  $seqObj->display_id;
        my $fullRefseqSeq = $self->{refseq}{$refseqID}{seq};

        ###################################################
        #Step1: Process Refseq alignment
        ###################################################
        my $isReferenceSeq  =  ($refseqID =~ m/^ref\|/);
        my $isGood = $self->{refseq}{$refseqID}{safe};
        if($isReferenceSeq && $isGood)
        {
            delete $self->{refseq}{$refseqID}{ntMSALoc};
            my $fullAlignment = $seqObj->seq;
            my $offset = 0;
            while($fullAlignment =~ m/(?<continuousAA>(?<codon>[^-]--)+)/g)
            {  #will only match X-- patterns so wont work for ----
                my $start           =  $-[0];  #zero indexed amino acid position
                my $spacedSeq       =  $+{continuousAA};
                (my $continuousSeq  =  $spacedSeq) =~ s/-//g;

                my $index           =  index(substr($fullRefseqSeq, $offset), $continuousSeq);  #incomplete alignment with the reference sequernce

                #within the single stretch process each AA
                my $aa              =  0;
                while($spacedSeq =~ m/[^-]/g)
                { #move stepwise from 1 aa to the next aa
                    $aa++;
                    my $ntLoc = $start + $-[0] + 1; #plus one cause ltr we will do substr
                    #msa_nt_LOCATION :: aaLOCATION
                    $self->{refseq}{$refseqID}{ntMSALoc}{$ntLoc} = $aa + $offset + $index;
                }
                $offset += length $continuousSeq;
            }

            ###################################################
            #Step2: NT Contigs
            ###################################################
            while($seqObj = $in->next_seq){
                my $contigID = $seqObj->display_id;
                my $isComplement = $self->{contigs}{$contigID}{parentREF} eq $refseqID;
                if($isComplement)
                {
                    my $aaXaa = 0;  #new sequence
                    my $contigSequence = $seqObj->seq;
                    while($contigSequence =~ m/(?<codon>[^-]{3})/g)
                    {
                        #say "Aln: $alignmentFile\nQuery: $contigID\nSubject: $refseqID\ncodon: $+{codon}";
                        ##position in contigXrefseq msa
                        my $msaLoc  =  $-[0] + 1;
                        #say "msa: ", substr($contigSequence, 0, $msaLoc+2);

                        ##position in aa
                        #Match to an amino acid may or may not exists
                        #Protein D--  ---  --- G--
                        #Contig  ACC (gtg) gcc GCC
                        my $aa;
                        if(exists $self->{refseq}{$refseqID}{ntMSALoc}{$msaLoc})
                        {
                            $aa = $self->{refseq}{$refseqID}{ntMSALoc}{$msaLoc};
                            #say "aa position: ", $aa;
                            #my $tempAA = $aa - 1;
                            #say "the aa: ",join "", (split("",$self->{refseq}{$refseqID}{seq}))[0..$tempAA];

                            #position of aa in global msa (based on aa)
                            $aaXaa = $self->{refseq}{$refseqID}{map}{$aa};
                            #store codon in contig obj for that aaXaa position
                            $self->{contigs}{$contigID}{globalCoordinates}{$aaXaa} = $+{codon};
                        }else
                        {
                            #if there isnt any match of aa then increment the position by 1/1e6
                            $aaXaa+=1/1e6;
                            #say "doesnt exists: $contigID\t$aaXaa";
                            $self->{contigs}{$contigID}{globalCoordinates}{$aaXaa} = $+{codon};
                        }
                        $self->{pos}{$aaXaa}++;
                    }
                }else{
                    next
                }
            }
        }else{
            next
        }
    }
}


#Accessory methods
sub grepRefSeqID ($self, $header){
$header =~ m/(ref\|\S+?\|)/                                                                                     ;
return $1;
}

sub minmax($self){
    #Criteria for selecting refseq sequences
    #       1. length must be within 2 SDs of the mean length left and right (will no)
    my @len = map {$self->{refseq}{$_}{length} }keys $self->{refseq}->%*                                                                      ;
    $self->min(mean(@len) - 2 * stddev(@len));
    $self->max(mean(@len) + 2 * stddev(@len));
    print "#\tMIN ref length: ",$self->min,"\n";
    print "#\tMAX ref length: ",$self->max,"\n";
}


has pos => (
    is => 'rw',
    handles => {},
);

#Curated
has refseqFasta => (
    is => 'ro',
    required => 1,
);

has alignmentFiles =>(
is=>'ro',
handle => [],
#required=>1,
);

has outputPrefix =>(
is=>'ro',
required=>1,
);

#Refseq sequences & lengths
has refseq => (
    is =>'rw',
    handle => {},
);

#Contigs
has contigs => (
is=>'rw',
handle=>{},
);

#Min and max lengths of refseq sequences
has min=>(
is=>'rw',
default=>0,
);

has max=>(
is=>'rw',
default=>0,
) ;
1;
