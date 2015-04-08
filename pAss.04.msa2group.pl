#!/usr/bin/env perl


use Modern::Perl '2014';
use autodie;
use Statistics::R;

die "usage: $0 megan.aln out.prefix\n" unless $#ARGV == 1;
my $in = shift;
my $prefix = shift;

my $kmer = 25;
my $step = 5;
my $too_much_gap = 10;
my $min_reads = 10;
my $skip = 0.1;
#my $gap_skip = int($kmer * $skip);
my $gap_skip = 4;

open my $IN, "<", $in;
my %aln = <$IN>;
my @seq = values %aln;
my $total = scalar @seq;
if($total < $min_reads)
{
	print "too few reads to process: $total\n";
	exit;
}
my $gaps;
my $nongaps;
foreach(@seq)
{
	chomp;
	s/^(-+)/'.' x length($1)/e;
	s/(-+)$/'.' x length($1)/e;
	$gaps += s/-/-/g;
	$nongaps += s/([ATGCNatgcn])/$1/g;
}
my $gap_ratio = $gaps / $nongaps;
if($gap_ratio >= $too_much_gap)
{
	print "too much gaps: $gaps gaps vs $nongaps non-gaps\n";
	exit;
}

open RAW, ">$prefix.csv";
say RAW "position,nseq,ngroup,total";
open DIS, ">$prefix.dist.csv";
say DIS "position,nseq,gsize";

my $aln_size = length($seq[0]);

for my $i(0..($aln_size - $kmer)/$step) #the number of windows
{
	my $nseq;
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

__END__
system(qq(echo 'd <- read.csv("$prefix.csv"); pdf("$prefix.2.pdf", w = 20, h = 8, pointsize=8); par(mfrow=c(2,1), mar=c(4,4,1,3)); with(d, plot(position, ngroup)); with(d, plot(position, nseq))' | R --vanilla));

