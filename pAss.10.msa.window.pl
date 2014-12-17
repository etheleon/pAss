#!/usr/bin/env perl

use Modern::Perl '2014';
use autodie;
use Bio::SeqIO;
use List::Util qw(sum max);
use Array::Utils qw(:all);
use Data::Alias;
use Method::Signatures;
use Set::IntervalTree;

die "usage: $0 min_leng[integer] msa.file\n" unless $#ARGV == 1;
##################################################
#+------------------------------------------------
#Init
#+------------------------------------------------
##################################################

my ($minLength, $inputFile) = @ARGV;    #cutRegionSize, input.msa
my $tree = Set::IntervalTree->new;
my $msalength;  #length of msa
my %seq;        #hash ref storing sequence information
my %nnucl;      #depth at postion i of MSA
my @starting;
$|++;
##################################################
#+------------------------------------------------
#Main
#+------------------------------------------------
##################################################
my ($KO) = $inputFile =~ m/(K\d{5})/;
#Step1: Parse MSA
my $input = Bio::SeqIO->new(-file => "$inputFile", -format=>"fasta");
while (my $seqObj= $input->next_seq)
{
    parseMSA($seqObj->display_id, $seqObj->seq)
    #stores - msa length
    #       - sequence information
    #       - build intervals
}

#Step2: Find 10bp window with MAX diversity
my $max_seq = max map {windowDepth($_, 10)} (9..$msalength);
my $nseq = scalar keys %seq;    #dunno what is this for
#Step3
#Calculate window statistics
say join "\t", qw(ko msaTotSeq start end maxSeq.10bp seqInSameWindow newCut lens);
slidingWindow($_) for (0..($msalength-$minLength));

##################################################
#+------------------------------------------------
#Functions
#+------------------------------------------------
##################################################

func parseMSA ($header,$seq)
{
    my $count = ($seq=~s/([^\-])/$1/g); #length without gaps
    buildInterval($header, $seq) unless $count < $minLength;
    $msalength = eval(length($seq) -1) unless defined $msalength;
}

func buildInterval ($header, $sequence)
{
    if($sequence =~ m/^ (?<frontGap>-*) (?<body>.+?)    (?<endGap>-*)$/x)
    {
        my $start = length $+{frontGap};
        my $end   = $start + length $+{body};

        $seq{$header} = $sequence;

        #Objective2: Build intervals
        $tree->insert($header, $start, $end);
        $nnucl{$_}++ for ($start..$end);
    }
}

func windowDepth ($position, $windowSize)
{
    my @winSize = (0..eval($windowSize-1));
    my $totsize = sum map { $nnucl{$position - $_} } @winSize;
    return eval($totsize /= 10);
}

#stores sequences found in the window and announces the max size
func slidingWindow ($startingPosition)
{
    @starting = ();
    @starting = @{$tree->fetch($startingPosition, $startingPosition+1)};
    #Stores headers of sequences in the starting position
    tryCutSize($minLength, $startingPosition);
}

func tryCutSize($cut, $startingPosition)
{
    my $endingPosition = $startingPosition + $cut;
    if ($endingPosition < $msalength)
    {
        #optimise
        my @ending   = @{$tree->fetch($endingPosition, $endingPosition+1)};
        my @isect    = intersect @starting, @ending;
        #optimise
        my $ns       = scalar @isect;
        if($ns)
        {
            my @lens = map { countNoGap($_, $startingPosition, eval($cut)) } @isect;
            #the averaged number of nucleotides (excluding gaps) across sequences captured in this window
            my $lens = sprintf "%.02f", (sum(@lens)) / (scalar @lens);
            if ($lens >= $minLength)
            {
                my @row = (
                $KO,                # the ko
                $nseq,              # the number of sequences in that MSA
                $startingPosition,  # the start of window
                $endingPosition,    # the end of window
                $max_seq,           # the num of sequences captured within window
                $ns,                # the number of sequences spanning window
                $cut,
                $lens);             # the average length of the spanning sequences excluding gaps
                say join "\t", @row;
            }else
            {
                $cut += int max(1, $minLength - $lens + 1);
                #add to cut window either 1 or num.gaps + 1
                tryCutSize($cut, $startingPosition)
            }
        }else{return}
    }else{return}
}

func countNoGap($header, $starting, $length)
{
    my $sequence = $seq{$header};
    my $this = substr($sequence, $starting, $length);
    my $nogap = ($this =~ s/[^-]//g);
    return $nogap;
}
