#!/usr/bin/env perl
use strict;
use List::Util qw(max);

my $MIN_OVERLAP = 40;
my $MIN_IDENTITY = 0.95;

my $prefix = shift;

my $out = shift;
my $log = shift;

die "usage: $0 prefix OUT include.list alignment-file1 alignment-file2 ...\n" unless @ARGV;
die "usage: $0 prefix OUT include.list alignment-file1 alignment-file2 ...\n" if $out =~ m/\/alignment-/;
die "usage: $0 prefix OUT include.list alignment-file1 alignment-file2 ...\n" if $log =~ m/\/alignment-/;

die unless @ARGV;
open OUT, ">$out" or die;
open LOG, ">$log" or die;

my @files = @ARGV;

my %seq;
my %overlap;
my %how;
my %best;
my %info;
my %included;

for my $file(@files)
{
#	print STDERR "$file\n";
    open IN, $file or die $!;
    my $temp = <IN>;
    $temp = <IN>;

    my %left;
    my %right;
    while(my $head = <IN>)
    {
        chomp $head;
        $head =~ s/ .+?(\(rev\))?$/$1/;
        $head =~ s/\(rev\)/-rc/;
        $head =~ s/>//;
        my $seq = <IN>;
        chomp $seq;
        $seq =~ uc $seq;

        my $left = $1 if $seq =~ m/^(-+)/;
        $left = 0 + (length $left);

        my $right = $1 if $seq =~ m/(-+)$/;
        $right = length($seq) - length($right) - 1;

        $seq{$file}->{$head} = $seq;
        $left{$head} = $left;
        $right{$head} = $right;
    }

# order MUST be from left to right, otherwise substr($seq, $lefty, $overlap) won't work
    my @reads = sort {$left{$a} - $right{$a}/1e6 <=> $left{$b} - $right{$b}/1e6} keys %left;
    for my $i(0..$#reads)
    {
        my $x = $reads[$i];
        my $left = $left{$x};
        my $right = $right{$x};
        my $seq = $seq{$file}->{$x};
#		print "\n\n$left\t$right\t$x\n";
        for my $j(($i+1)..$#reads)
        {
            my $y = $reads[$j];
            my $lefty = $left{$y};
            my $righty = $right{$y};

            my $overlap = ($right - $left + 1) - max($lefty - $left, 0) - max($right - $righty, 0);
            next unless $overlap >= $MIN_OVERLAP;

            my $seqy = $seq{$file}->{$y};

            my $xx = substr($seq,  $lefty, $overlap);
            my $yy = substr($seqy, $lefty, $overlap);
            my $match = match($xx, $yy);

            if($match / $overlap >= $MIN_IDENTITY)
            {
#				print "\t$lefty\t$righty\t$y\n";
#				print "\toverlap: $overlap\n";
#				print "\tmatch: $match\n";
#				print "\tpercent: ", (100 * $match / $overlap), "\n";
                if($right >= $righty)
                {
                    unless($included{$y})
                    {
                        print LOG "$y\n";
                        $included{$y} = 1;
                    }
                    next;
                }

                if($match > $overlap{$x}->{$y})
                {
                    $overlap{$x}->{$y} = $match;
                    $overlap{$y}->{$x} = $match;

                    my $seqx = $seq{$file}->{$x};
                    my $seqy = $seq{$file}->{$y};
                    my $seq0 = substr($seqx, 0, $lefty);
                    my $seq2 = substr($seqy, $lefty+$overlap);

                    my $seq1;
                    my $seqxx = $seqx;
                    my $seqyy = $seqy;
                    $seqxx =~ s/-//g;
                    $seqyy =~ s/-//g;
                    my $lenx = length $seqxx;
                    my $leny = length $seqyy;
                    if($lenx >= $leny)
                    {
                        $seq1 = substr($seqx, $lefty, $overlap);
                    }else
                    {
                        $seq1 = substr($seqy, $lefty, $overlap);
                    }
                    my $seq = $seq0.$seq1.$seq2;
                    $seq =~ s/-//g;
                    my $len = length $seq;

                    if($len > length($how{$x}->{$y}) && $len > $best{$y})
                    {
                        $how{$x}->{$y} = $seq;
                        $best{$y} = $len;

                        $info{$x}->{$y} = "x=$lenx:$left-$right y=$leny:$lefty-$righty overlap=$overlap match=$match";
                    }
                }
            }
        }
    }
}

my %done;
my $i= 0;
for my $x(keys %how)
{
    for my $y(sort {length $how{$x}->{$b} <=> length $how{$x}->{$a} } keys %{$how{$x}})
    {
        next if $done{$y}->{$x};

        my $len = length $how{$x}->{$y};
        next if $len < $best{$y};

        $i++;
        print OUT ">$prefix.ass$i $x+$y length=$len $info{$x}->{$y}\n$how{$x}->{$y}\n";

        $done{$x}->{$y} = 1;
        last;
    }
}

sub match
{
    my $xx = shift;
    my $yy = shift;
#	print "$xx\n";
#	print "$yy\n";
    my $mask = "\xFF" x length($xx);
    $xx ^= $yy;
    $xx &= $mask;
#	$xx =~ s/[^\x00]/1/g;
#	$xx =~ s/\x00/0/g;
    my $match = ($xx =~ s/\x00//g);
    return $match;
}

__END__


$overlap_size = ($gene_right - $gene_left + 1) - max($left - $gene_left, 0) - max($gene_right - $right, 0);
