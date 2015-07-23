#!/usr/bin/env perl

use Modern::Perl '2015';
use autodie;
use Statistics::R;
use experimental qw/signatures/;
use Getopt::Lucid qw( :all );

my @specs = (
    Param("msa|m"),
    Param("outputPrefix|o"),
    Param("kmer|k")->default(25),
    Param("step|s")->default(5),
    Param("gap|g")->default(10),
    Param("minReads|r")->default(10),
    Param("skip|p")->default(0.1),
    Param("gapSkip|a")->default(4),
    Switch("help|h")->default(0)
);
my $opt = Getopt::Lucid->getopt( \@specs );
pod2usage(-verbose=>2) if $opt->get_help;

my $in     = $opt->get_msa;
my $prefix = $opt->get_outputPrefix;

#Arbituary thresholds
my $kmer         = $opt->get_kmer;
my $step         = $opt->get_step;
my $too_much_gap = $opt->get_gap;
my $min_reads    = $opt->get_minReads;
my $skip         = $opt->get_skip;
my $gap_skip     = $opt->get_gapSkip;

open my $IN, "<", $in;
my %aln   = <$IN>; #key  = ID; value = sequence
my @seq   = values %aln;
my $total = scalar @seq;

# Check min reads
if($total < $min_reads)
{
	say STDERR "too few reads to process: $total\n EXITING...";
	exit;
}

#Check Gaps
my $gaps;
my $nongaps;
foreach(@seq)
{
    chomp;
    s/^(-+)/'.' x length($1)/e; #trim head
    s/(-+)$/'.' x length($1)/e; #trim tail
    $gaps += s/-/-/g;
    $nongaps += s/([ATGCNatgcn])/$1/g;
}

my $gap_ratio = $gaps / $nongaps;
if($gap_ratio >= $too_much_gap)
{
	say "too much gaps: $gaps gaps vs $nongaps non-gaps\n EXITING";
	exit;
}

open RAW, ">$prefix.csv";
say RAW "position,nseq,ngroup,total";

open DIS, ">$prefix.dist.csv";
say DIS "position,nseq,gsize";

my $aln_size = length($seq[0]);
for my $i(0..($aln_size - $kmer)/$step) #the number of windows
{
	my $nseq = 0;
	my %group;
	my $pos = $i * $step;
	for my $seq(@seq)
	{
        my $mer = substr($seq, $pos, $kmer);
        next if $mer =~ m/\./;
        next if $mer =~ m/-{$gap_skip,}/;
        $group{$mer}++;
        $nseq++;
	}
	my $ngroup = keys %group;
	if($nseq >= $min_reads)
	{
		print RAW "$pos,$nseq,$ngroup,$total\n";
		for my $size(values(%group))
		{
			my $prop = $size / $nseq;
			print DIS "$pos,$nseq,$prop\n";
		}
	}
}

my $R = Statistics::R->new();
my ($lines) = `wc -l $prefix.csv` =~ m/^(\d+)/;
if($lines > 1)
{
    $R->run(qq`attach(read.csv("$prefix.csv"))`);
    $R->run(qq`pdf("$prefix.1.pdf");plot(nseq, ngroup);abline(0,1);dev.off()`);
    $R->run(qq`pdf("$prefix.2.pdf", w = 10); par(mfrow=c(2,1), mar=c(4,4,1,3)); plot(position, ngroup); plot(position, nseq)`);
}

=pod

=head1 NAME
    Some Name

=head1 SYPNOSIS

    Some SYPNOSIS

=head1 DESCRIPTION

    The reason behind the options kmer, step, gap, minReads, skip, gapSkip are entirely arbituary. Original author Xie Chao has not documented why this was so

=head1 OPTIONS

=over 4

=item --msa -m

    nucleotide alignment from pAss03

=item --outputPrefix -o

    prefix (excluding the file extension) to the output files

=item --kmer -k

    kmer size; default 25

=item --step -s

    step size; default 5

=item --gap -g

   max allowed gaps; default 10

=item --minReads -r

    minimum allowed reads; default 10

=item --skip -p

    skip threadhold; dont really know what this means; default 0.1

=item --gapSkip -a

    default 4; the number of gaps to skip

=back
=cut

__END__
OLD code
system(qq(echo 'd <- read.csv("$prefix.csv"); pdf("$prefix.2.pdf", w = 20, h = 8, pointsize=8); par(mfrow=c(2,1), mar=c(4,4,1,3)); with(d, plot(position, ngroup)); with(d, plot(position, nseq))' | R --vanilla));
