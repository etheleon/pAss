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

# ABSTRACT: Gene centric diversity indexing

=pod

=head1 NAME

pAss library

=head1 SYNOPSIS

library for generating prot MSA guided MSA of nt contigs

=head1 Initialise

C<my $alignment = PASS::Alignment->new(>
C<    refseqFasta     =>  "$FindBin::Bin/data/refSeqProtDB/ko\:K00001", #reference sequences>
C<    alignmentFiles  =>  [glob("$FindBin::Bin/data/pAss01/K00001/alignment-*")],>
C<    outputPrefix    =>  "$FindBin::Bin/data/pAss",>
C<)>

=head1 Methods

=cut

=pod

=over

=item $obj->storeRefseq()

Stores information about the protein reference sequences.
Calculates the min and max length allowed for the sequence to be included in the MSA

=back

=cut

sub storeRefseq($self){
    my $dbin = Bio::SeqIO->new(
        -file=>$self->refseqFasta,
        -format=>"fasta");
    while(my $seqObj = $dbin->next_seq)
    {
        my $refseqID = $self->grepRefSeqID($seqObj->display_id);
        $self->{refseq}{$refseqID} = {seq => $seqObj->seq, length => $seqObj->length}
    };
    $self->minmax;
	for my $refseqID (keys $self->{refseq}->%*){
		my $isCorrectLength = $self->{refseq}{$refseqID}{length} > $self->min &&  $self->{refseq}{$refseqID}{length} < $self->max;
        $self->{refseq}{$refseqID}{safe} = $isCorrectLength ? 1 : 0;
	}
}

=pod

=over

=item $obj->assignContig2ref( %options )

=head3 Determines which contig be mapped to which protein reference sequences

In some cases, no contigs are considered for a particular KO because:

=over 3

=item 1. length of sequences

Entirely KO is left out cause their reference sequences are all the same length;
SD is low and the alignment length is not the length (this doesnt make sense)

=item 2. # of observations

Too little reference sequences to build a MSA.

=back

For the same sequence, we observe it being put into different files by MEGAN. eg. with different alignment...
NOTE: this is already taken care of by assignContig2ref so no problem

data/pAss01/K00001/alignment-Leptospira_interrogans-ref_NP_711797.1__alcohol_dehydrogenase__Leptospira_interrogans_serovar_Lai_str.-01029.fasta
>ref|NP_711797.1| alcohol dehydrogenase [Leptospira interrogans serovar Lai str.
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------A--K--A--L--S--V--M--M--C--N--E--G--S--A--S--Q--Y--Q--G--W--D--L--V--P--N--P--G--I
--Y--K--I--G--I--P--T--V--A--G--S--G--A--E--A--S--R--T--A--V--L-----------------------------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------V--F--L--S--D--G--D--D--D--M--L--L--T--A--S--Y--M--G--G--V--S--I--V--N--S--E--V--G--V--C--H--A--L--S--Y--G--L--S--L--E--L-
-G--Y--R--H--G--F--A--


data/pAss01/K00001/alignment-cellular_organisms-ref_NP_711797.1__alcohol_dehydrogenase__Leptospira_interrogans_serovar_Lai_str.-01553.fasta
>ref|NP_711797.1| alcohol dehydrogenase [Leptospira interrogans serovar Lai str.
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------G--W--D--L--V--P--N--P--G--I
--Y--K--I--G--I--P--T--V--A--G--S--G--A--E--A--S--R--T--A--V--L--M--G--K--E--R--K--F--G--I--N--S--D--H--S--M--F--D--A--I--I--L--D--S--S--L--I--K--N--V--P--I--A--Q--R--F--Y--S--G--M--D--C--Y--I--H--C--
V--E--S--L--Q--G--T--M--I--N--E--L--A--K--G--N--A--S--K--A--L--E--L--C--E--K--V--F--L--S--D--G--D--D--D--M--L--L--T--A--S--Y--M--G--G--V--S--I--V--N--S--E--V--G--V--C--H--A--L--S--Y--G--L--S--L--E--L-
-G--Y--R--H--G--F--A--N--C--V--A--F--N--V--L--D--E--Y--Y--G--P--W--V--D--R--F--R--E--M--L--K--I--H--

=back

=cut

sub assignContig2ref($self){
    for my $alignmentFile ($self->alignmentFiles->@*)
    {
        my $in = Bio::SeqIO->new(-file=>$alignmentFile,-format=>'fasta');
			#INJEST PROTEIN SEQUENCE
			my $refObj = $in->next_seq;
			my $refseqID = $self->grepRefSeqID($refObj->display_id);
			next unless $self->{refseq}{$refseqID}{safe};

		#PARSE NT CONTIG
        while(my $seqObj = $in->next_seq)
		{
            my $contigID  =  $seqObj->display_id;
            my $contigDESC = $seqObj->desc;
            my $direction = $contigDESC =~ m/\(rev\)/ ? '(rev)' : '(ntRev)';
            my $sequence = $seqObj->seq;
            my $gap = ($sequence =~ s/-//g);
            my $length = length $sequence;
            $length -= $gap / 1e6; #this parameter represents the length of the alignment - the number of gaps,
            #say "$alignmentFile\n$refseqID\t$length" if $contigID eq 'contig00049';
            my $contigDetails = {
                seq  	  =>  $seqObj->seq,
                len 	  =>   $length,
                parentREF => $refseqID,
                direction => $direction
            };
            if(!exists $self->{contigs}{$contigID}){
                $self->{contigs}{$contigID} = $contigDetails;
            }else{
				#if alignment is longer, then use this alignment
                my $isBetterMatch = $contigDetails->{len} > $self->{contigs}{$contigID}{len};
                if ($isBetterMatch){
                    #say "Found better match";
                    $self->{contigs}{$contigID} = $contigDetails
                }
			}
		}
	}
}

=pod

=over

=item $obj->runMuscle( %options )

Runs muscle on the stored reference protein sequence.

Whether to load a fixed set of aligned MSAs or run a MUSCLE from the beginning

=back

=cut

sub runMuscle($self, $preran=0){
    #muscle is heuristic, each run will give you a different output
		my $out = $self->outputPrefix;
	unless($preran){
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
	}

	my $inputFile = $preran ? $preran : "$out.temp.ref.msa";
    my $MUSCLEinput=Bio::SeqIO->new(-file=>$inputFile, -format=>"fasta");
    while(my $seqObj = $MUSCLEinput->next_seq)
	{
        my $seq = $seqObj->seq;
        my $refseqID = $seqObj->display_id;
        my $aaIndex = 1;
        while($seq =~ m/[^-]/g)
        {
            #$-[0] is zero index
            my $msaLOC = $-[0] + 1;
            $self->{refseq}{$refseqID}{map}{$aaIndex} = $msaLOC;
            $aaIndex++;
        }
    }
}

=pod

=over

=item $obj->buildContigMSA()

Two step process,
Stores the location of each amino acid from the blastX alignment in the alignment object.

Process

	1. loops thru each contig-reference protein alignment (rmbr there might be duplicates, ie above)
		CHECK: refseq prot meets the size requirement and alignment.
	2.

=back

=cut

sub buildContigMSA($self){
    #gives a reference keypair of locations key:position value:aa
    foreach my $alignmentFile ($self->alignmentFiles->@*){
        my $in=Bio::SeqIO->new(-file => $alignmentFile, -format => "fasta");
        my $seqObj = $in->next_seq;
        my $refseqID = $self->grepRefSeqID($seqObj->display_id);

        my $isReferenceSeq  =  ($refseqID =~ m/^ref\|/);
        my $isGood = $self->{refseq}{$refseqID}{safe};
        if($isReferenceSeq && $isGood)
        {
		#i think i have to rewrite this portion OMG!!!!
        ###################################################
        #Step1: Process Refseq alignment
        ###################################################
			my $fullRefseqSeq = $self->{refseq}{$refseqID}{seq};
            delete $self->{refseq}{$refseqID}{ntMSALoc}; #the ntMSALoc is only temporary for that file; the really impt one is for
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
            #Step2: NT Contigs: its not just NT sequences from a single alignment file handed down from
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
                            #if there isnt any match of aa ie. the codons match to a gap then increment the position by 1/1e6
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

#File path pointing to file storing protein sequences of tall reference sequences
has refseqFasta => (
    is => 'ro',
    required => 1,
);

#MSA alignment files created by MEGAN5 to of contigs against reference sequences
has alignmentFiles =>(
is=>'ro',
handle => [],
#required=>1,
);

#the output directory and prefix
has outputPrefix =>(
is=>'ro',
required=>1,
default=>"out"
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
