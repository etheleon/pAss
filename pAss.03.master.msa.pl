#!/usr/bin/env perl

use v5.20;
use experimental 'signatures';
use autodie;
use Bio::SeqIO;
use msa::alignment;

#Pragma
use strict;
use warnings;

use Statistics::Basic qw(mean stddev);
use Set::Intersection;

my ($dbfasta, $msadir, $out) = @ARGV;
die "usage: $0 db.fasta msa.dir out.msa\n" unless $#ARGV == 2 && -f $dbfasta && -d $msadir;
$msadir =~ s|//|/|g; $msadir =~ s|/$||; #deal with weird path names

my %seq;        #stores reference sequences
my %good_hit;   #stores contigs which are good matches to the reference sequences
my $temp = "$out".rand().time();

##################################################
say STDERR "#reading db file";
##################################################
my $in = Bio::SeqIO->new(-file=>$dbfasta, -format=>"fasta");
while(my $seqObj = $in->next_seq){
    my $seq = $seqObj->seq();
    $seq =~ s/\s+//g;
    my $id = $seqObj->display_id();
    $id = $1 if $id =~ m/(ref\|.+?\|)/;
    $seq{$id}{seq}     =  uc $seq;
    $seq{$id}{length}  =  length $seq;}

my @len = map $seq{$id}{len} keys %seq;
my $len = { min=>mean(@len) - 2 * stddev(@len),
            max=>mean(@len) + 2 * stddev(@len)};

##################################################
say STDERR "#reading list of good references";
##################################################
selectRefseqs($_) for glob "$msadir/alignment-*";

##################################################
say STDERR "#preparing good_ref sequences for MSA";
##################################################
open my $out, ">", "$out.temp.ref.faa";
my %good_ref;
for my $ref (keys %good_hit){
    next if $good_ref{$good_hit{$ref}{id}};
    print $out ">$good_hit{$ref}{id}\n$seq{$good_hit{$ref}{id}}\n";
    $good_ref{$ref} = 1;}

##################################################
say STDERR "#Running MSA with Muslce";
##################################################
system "muscle -in $out.temp.ref.faa -out $out.temp.ref.msa";

##################################################
say STDERR "#Mapping ref_seq to MSA coordinates";
##################################################
my %map;
my $in=Bio::SeqIO->new(-file=>"$out.temp.ref.msa", -format=>"fasta");
while(my $seqObj = $in->next_seq){
        my $id = $seqObj->display_id;
        my $seq = $seqObj->seq;
        $seq =~ s/\s+//g;
        while($seq =~ m/[^-]/g){
            state $aaIndex++;
            my $subStrLoc= $-[0] + 1;
            $map{$id}{$aaIndex} = $subStrLoc;}}}

##################################################
say STDERR "mapping raw DNA sequences to MSA coordinates";
##################################################

my %hit;
my %msa2ref;
my $ref;
my %pos;

##################################################
print STDERR "$msaf\n";
##################################################
loopMSAF($_) for glob "$msadir/alignment-*";

open OUT, ">","$out.msa";
my $out=Bio::SeqIO->new(-file=>">out.msa",-format=>"fasta")
for my $id (keys %hit)
{
    say OUT ">$id";
    for my $pos(sort {$a <=> $b} keys %pos)
    {
        my $codon = $hit{$id}->{$pos};
        $codon = '---' unless $codon;
        print OUT $codon;
    }
    print OUT "\n";
}
close OUT;

sub selectRefseqs ($alignment)                                  {
    say STDERR "$msaf"                                          ;
    my $in=Bio::SeqIO->new(-file=>$alignment, -format=>'fasta') ;
#PROTEIN SEQUENCE
    my $refObj = $in->next_seq                                  ;
    my $head   = $refObj->display_id                              ;
    if($head =~ m/^>(\S+)/)                                     {
        #skip if the sequence is too short
        my $ref = $1                                            ;
        return if $seq{$ref}{length} < $len{min}                 ;
        return if $seq{$ref}{length} > $len{max}                 ;}
#Parse Contig MSA (DNA)
    while(my $seqObj = $in->next_seq)                           {
        my $h = $seqObj->display_id
        $h = $1 if $h =~ m/>(\S+)/                              ;
        my $gap = ($seq =~ s/-//g)                              ;
        my $len = $seqObj->length                               ;
        $len -= $gap/1e6                                        ;
        if($len > $hit_len{$h})                                 {
            $good_hit{$h}{id} = $ref                            ;
            $good_hit{$h}{length} = $len                        ;}}}

sub loopMSAF ($alignmentFile){
my $in=Bio::SeqIO(-file=>$msaf,-format=>"fasta");
    while(my $seqObj = $in->next_seq){
        my $id = $seqObj->display_id;
        my $seq = $seqObj->seq;
        $seq = uc $seq;
        $seq =~ s/\s+//g;
        my $alignment = msa::alignment->new(sequence=>$seq, id=>$id):
        if($alignment->ifRefseq)
        {
            while($seq =~ m/(([^-]--)+)/g)
            {
                my $start = $-[0];
                my $fullseq = $1;
#			print "$stripseq\n";
                my $index = index(substr($seq{$id}, $last_end), $stripseq);
#			print "$index\n";
                my $i = 0;
                while($fullseq =~ m/[^-]/g)
                {
                    $i++;
                    my $j = $start + $-[0] + 1;
                    $msa2ref{$j} = $i + $last_end + $index;
#						print "\t\t", substr($seq{$id}, $msa2ref{$j}-1, 1), "=", substr($seq, $j-1, 1), "\n";
                }
                $last_end += length $stripseq;
            }
        }else
        {
            next unless $good_hit{$id} eq $ref;

            while($seq =~ m/([^-]{3})/g)
            {
                my $codon = $1;
                my $j = $-[0] + 1;
                my $i = $msa2ref{$j};
                my $this = $map{$ref}{$i};
                if($this)
                {
                    $pos = $this;
                }else
                {
                    $pos += 1/1e6;
                }
                $pos{$pos}++;
                $hit{$id}->{$pos} = $codon;
            }
        }
    }
}


}

package msa::alignment;

use Moo;
use v5.20;
use experimental 'signatures';
use strictures 2;
use namespace::clean;

sub isGoodSeq ($self){
    my $id = $self->id;
    my ($isRefseq) = $id =~ m/^ref\|/ ? 1 : 0;
    return exists $good_ref{$ref} ? 1 : 0;
}

has stripSeq =>(
    is      => 'ro',
    default => sub($self){
        my $sequence = $self->sequence;
        return $sequence =~ s/-//g;
    }
)

has sequence => (
    is => 'ro',
    isa => sub {die "No spaces allow" unless $_[0] =~ !/\s*/},
    required=>1,
);

has position => (
    is  => 'rw',
    default=> 0,
);

has last_end => (
    is =>'rw',
    default=>0,
);

has id => (
    is  => 'ro',
    required=>1,
);



1;

if($id =~ m/^ref\|/)
    {
        $ref = $id;
        next unless $good_ref{$ref};

        my $last_end = 0;
        while($seq =~ m/(([^-]--)+)/g)
        {
            my $start = $-[0];

            my $fullseq = $1;
            my $stripseq = $1;
            $stripseq =~ s/-//g;
#					print "$stripseq\n";

            my $index = index(substr($seq{$id}, $last_end), $stripseq);
#					print "$index\n";

            my $i = 0;
            while($fullseq =~ m/[^-]/g)
            {
                $i++;
                my $j = $start + $-[0] + 1;
                $msa2ref{$j} = $i + $last_end + $index;
#						print "\t\t", substr($seq{$id}, $msa2ref{$j}-1, 1), "=", substr($seq, $j-1, 1), "\n";
            }


            $last_end += length $stripseq;
        }
    }else
    {
        next unless $good_hit{$id} eq $ref;

        while($seq =~ m/([^-]{3})/g)
        {
            my $codon = $1;
            my $j = $-[0] + 1;
            my $i = $msa2ref{$j};
            my $this = $map{$ref}{$i};
            if($this)
            {
                $pos = $this;
            }else
            {
                $pos += 1/1e6;
            }
            $pos{$pos}++;
            $hit{$id}->{$pos} = $codon;
        }
    }
}
    }
}


