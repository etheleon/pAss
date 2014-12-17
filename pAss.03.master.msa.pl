#!/usr/bin/env perl

use strict;
use lib "/export2/home/uesu/perl5/lib/perl5";
use Statistics::Basic qw(mean stddev);
use Set::Intersection;

my ($dbfasta, $msadir, $out) = @ARGV;
die "usage: $0 db.fasta msa.dir out.msa\n" unless $#ARGV == 2 && -f $dbfasta && -d $msadir;
$msadir =~ s|//|/|g;
$msadir =~ s|/$||;

my $temp = "$out".rand().time();

print STDERR "reading db file\n";
my %seq;
my %len;
open(IN, $dbfasta) or die;
my $term = $/;
$/ = "\n>";
while(<IN>)
{
    s/^\s+|\s+$//g;
    s/^\s*>|>\s*$//g;
    s/^\s+|\s+$//g;
    if(m/^(\S+).*?\n(.+)/s)
    {
        my $id = $1;
        my $seq = $2;

        $id = $1 if $id =~ m/(ref\|(.+?)\|)/;
        $seq =~ s/\s+//g;
        $seq{$id} = uc $seq;
        $len{$id} = length $seq;
    }
}
$/ = $term;

my @len = values %len;
my $len_mean = mean(@len);
my $len_sd = stddev(@len);
my $len_min = $len_mean - 2 * $len_sd;
my $len_max = $len_mean + 2 * $len_sd;

print STDERR "reading list of good references\n";
my %good_hit;
my %hit_len;
my $ref;
#my $xx;
for my $msaf(glob("$msadir/alignment-*"))
{
#	$xx++;
#	last if $xx > 30;
    print STDERR "$msaf\n";
    open(IN, $msaf) or die;
    my $head = <IN>;
    my $temp = <IN>;
    if($head =~ m/^>(\S+)/)
    {
        $ref = $1;
        next if $len{$1} < $len_min;
        next if $len{$1} > $len_max;
    }
    while(my $h = <IN>)
    {
        my $seq = <IN>;

        $h = $1 if $h =~ m/>(\S+)/;
        my $gap = ($seq =~ s/-//g);
        my $len = length $seq;
        $len -= $gap/1e6;
        if($len > $hit_len{$h})
        {
            $good_hit{$h} = $ref;
            $hit_len{$h} = $len;
        }
    }
}

print STDERR "preparing good_ref sequences for MSA\n";
open(OUT, ">$out.temp.ref.faa") or die;
my %good_ref;
for my $ref(values %good_hit)
{
    next if $good_ref{$ref};
    print OUT ">$ref\n$seq{$ref}\n";
    $good_ref{$ref} = 1;
}
print STDERR "MSA\n";
system("muscle -in $out.temp.ref.faa -out $out.temp.ref.msa");

print STDERR "mapping ref_seq to MSA coordinates\n";
my %map;
open(IN, "$out.temp.ref.msa") or die;
$/ = "\n>";
while(<IN>)
{
    s/^\s+|\s+$//g;
    s/^\s*>|>\s*$//g;
    s/^\s+|\s+$//g;
    if(m/^(\S+).*?\n(.+)/s)
    {
        my $id = $1;
        my $seq = $2;
        $seq =~ s/\s+//g;

        my $i = 0;
        while($seq =~ m/[^-]/g)
        {
            $i++;
            my $j = $-[0] + 1;
            $map{$id}->{$i} = $j;
        }
    }
}
$/ = $term;

print STDERR "mapping raw DNA sequences to MSA coordinates\n";
my %hit;
my %msa2ref;
my $ref;
my %pos;
#my $xx;
for my $msaf(glob("$msadir/alignment-*"))
{
#	$xx++;
#	last if $xx > 30;
    print STDERR "$msaf\n";
    open(IN, $msaf) or die;
    $/ = "\n>";
    while(<IN>)
    {
        s/^\s+|\s+$//g;
        s/^\s*>|>\s*$//g;
        s/^\s+|\s+$//g;
        if(m/^(\S+).*?\n(.+)/s)
        {
            my $id = $1;
            my $seq = $2;
            $seq = uc $seq;
            $seq =~ s/\s+//g;

            my $pos = 0;

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
                    my $this = $map{$ref}->{$i};
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

open OUT, ">$out.msa" or die;
for my $id(keys %hit)
{
    print OUT ">$id\n";
    for my $pos(sort {$a <=> $b} keys %pos)
    {
        my $codon = $hit{$id}->{$pos};
        $codon = '---' unless $codon;
        print OUT $codon;
    }
    print OUT "\n";
}
close OUT;

