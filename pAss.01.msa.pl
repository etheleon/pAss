#!/usr/bin/env perl
use strict;

die "usage: $0 query_file blast_file output.dir\n" unless $#ARGV == 2;

my $query = shift;
my $blast = shift;
my $out = shift;

die unless -f $query;
die unless -f $blast;

my $temp = rand().time();
my $CMD = join('', <DATA>);
$CMD =~ s/blastx.txt/$blast/;
$CMD =~ s/blastx.rma/$temp.rma/;
$CMD =~ s/query.fna/$query/;
$CMD =~ s/example/$out/;
open CMD, ">$temp" or die;
print CMD $CMD;
close CMD;
mkdir $out;


my $i;
my $signal = -1;
while($signal != 0)
{
    $i++;
    my $scr = int(30000 * rand());
    while(-e "/tmp/.X$scr-lock")
    {
        $scr = int(30000 * rand());
    }
    unlink "$out.lock";
    unlink "$out.log";
    unlink "$out.err";
    #$signal = system("xvfb-run -n $scr -f $out.lock -e $out.log ~/temp/megan/MEGAN +g -d -E < $temp");
    $signal = system("xvfb-run -n $scr -f $out.lock -e $out.log ~/local/megan/MEGAN -g -d -E < $temp");
    unless($signal == 0)
    {
        system("mv $out.log $out.err");
        system("echo $signal >> $out.err");
    }
    unlink "$out.lock";
    last if $i > 20;
}
unlink $temp;
unlink "$temp.rma";


__END__
import blastFile='blastx.txt' fastaFile='query.fna' meganFile='blastx.rma' maxMatches=25 minScore=50.0 maxExpected=1.0 topPercent=100.0 minSupport=1 minComplexity=0.30 useSeed=false useCOG=false useKegg=false paired=false useIdentityFilter=false textStoragePolicy=1 blastFormat=BlastX mapping='Taxonomy:BUILT_IN=true';
unCollapse nodes=all;
select nodes=all;
export what=alignment file='example/alignment-%c-%r-%n.fasta' data=Taxonomy classId=selected asConsensus=false useEachReadOnlyOnce=false useEachReferenceOnlyOnce=false refSeqOnly=false translateCDNA=false saveDiversityExtrapolation=false;
quit;
