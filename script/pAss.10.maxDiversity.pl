#!/usr/bin/env perl

use Modern::Perl '2015';
use autodie;
use Bio::SeqIO;
use experimental qw/signatures postderef/;
use List::Util qw(sum max);
use Array::Utils qw(:all);
use Set::IntervalTree;
use Pod::Usage;
use Getopt::Lucid qw/:all/;

#die "usage: $0 min_leng[integer] msa.file outputDir\n" unless $#ARGV == 2;
##################################################
#+------------------------------------------------
#Init
#+------------------------------------------------
##################################################

my @specs = (
    Param("minLength|m")->default(200),
    Param("inputFile|i"),
    Param("outputDir"),
    Switch("help|h")
);
my $opt = Getopt::Lucid->getopt( \@specs )->validate;
pod2usage(-verbose=>2) if $opt->get_help;

#my ($minLength, $inputFile, $outputDir) = @ARGV;    #cutRegionSize, input.msa
my ($minLength, $inputFile, $outputDir) = ($opt->get_minLength, $opt->get_inputFile, $opt->get_outputDir);

my $msalength;  #length of msa
my %seq;        #hash ref storing sequence information
my %nnucl;      #depth at postion i of MSA
my @starting;
my $tree    = Set::IntervalTree->new;
my ($KO)    = $inputFile =~ m/(K\d{5})/;

##################################################
#+------------------------------------------------
#Main
#+------------------------------------------------
##################################################

#Step1: Parse MSA
my $input = Bio::SeqIO->new(-file => "$inputFile",-format =>"fasta");

while (my $seqObj = $input->next_seq)
{
    parseMSA($seqObj->display_id, $seqObj->seq)
}

#Step2: Find 10bp window with MAX diversity
my $max_seq = max map {windowDepth($_, 10)} (9..$msalength);
my $nseq    = scalar keys %seq;

#Step3
mkdir $outputDir unless -d $outputDir;
my $outputFile = "$outputDir".'/'.$KO;
open my $output, ">", "$outputFile";
#Calculate window statistics
say $output join "\t", qw(ko msaTotSeq start end maxSeq.10bp seqInSameWindow newCut lens);
slidingWindow($_) for (0..($msalength-$minLength));

##################################################
#+------------------------------------------------
#Functions
#+------------------------------------------------
##################################################

sub parseMSA ($header,$seq)
{
    my $count   = ($seq=~s/([^\-])/$1/g); #length without gaps
    $msalength  = eval(length($seq) -1) unless defined $msalength;
    buildInterval($header, $seq) unless $count < $minLength;
}

sub buildInterval ($header, $sequence)
{
    if($sequence =~ m/^ (?<frontGap>-*) (?<body>.+?)    (?<endGap>-*)$/x)
    {
        my $start     = length $+{frontGap};
        my $end       = $start + length $+{body};
        say join "\t", $header, $start, $end;
        $seq{$header} = $sequence;

        #Objective2: Build intervals
        $tree->insert($header, $start, $end);
        $nnucl{$_}++ for ($start..$end);
    }
}

sub windowDepth ($position, $windowSize)
{
    my @winSize = (0..eval($windowSize-1));
    my $totsize = sum map { $nnucl{$position - $_} } @winSize;
    return eval($totsize /= 10);
}

#stores sequences found in the window and announces the max size
sub slidingWindow ($startingPosition)
{
    @starting = ();
    @starting = $tree->fetch($startingPosition, $startingPosition+1)->@*;
    #Stores headers of sequences in the starting position
    tryCutSize($minLength, $startingPosition);
}

sub tryCutSize($cut, $startingPosition)
{
    my $endingPosition = $startingPosition + $cut;
    if ($endingPosition < $msalength)
    {
        #optimise
        my @ending   = $tree->fetch($endingPosition, $endingPosition+1)->@*;
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
                say $output join "\t", @row;
            }else
            {
                $cut += int max(1, $minLength - $lens + 1);
                #add to cut window either 1 or num.gaps + 1
                tryCutSize($cut, $startingPosition)
            }
        }else{return}
    }else{return}
}

sub countNoGap($header, $starting, $length)
{
    my $sequence = $seq{$header};
    my $this     = substr($sequence, $starting, $length);
    my $nogap    = ($this =~ s/[^-]//g);
    return $nogap;
}
