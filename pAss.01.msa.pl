#!/usr/bin/env perl

use Modern::Perl '2015';
use autodie;
use experimental qw/signatures/;

die "usage: $0 query_file blast_file output.dir MEGAN\n" unless $#ARGV == 3;

my ($query, $blast, $out, $megan) = @ARGV;

die unless -f $query;
die unless -f $blast;
die unless -f $megan;

my $temp = rand().time();

prep($blast, $temp, $query, $out, $megan);
runMegan($src, $out, $temp, $stdFile, $megan, 1);

sub runMegan ($src, $out, $temp, $megan, $count){
    $count < 20 ? $count++ : return;

    unlink "$out.lock" if -e "$out.lock";
    unlink "$out.log" if -e "$out.log";

    my $scr = int(30000 * rand());

    #try new screen number if its taken;
    $scr = int(30000 * rand()) while -e "/tmp/.X$scr-lock";

    my $signal = `xvfb-run -n $scr -f $out.lock -e $out.log $megan -g -d -E  -c $temp`;

    #cant check for xvfb-run's own error [xc's version]
    $signal =~ m/Writing/sm ? eval{unlink "$out.lock", $temp, "$temp.rma": return} : runMegan($src, $out, $temp, $stdFile, $megan, $count);
}

sub prep ($blast, $temp, $query, $out, $megan){
    my $CMD = join '', <DATA>;
    $CMD =~ s/blastx.txt/$blast/;   #text file
    $CMD =~ s/blastx.rma/$temp.rma/;
    $CMD =~ s/query.fna/$query/;
    $CMD =~ s/example/$out/;

    open my $cmd, ">$temp" or die;
    print $cmd $CMD;
    close my $cmd;

    mkdir $out unless -d $out;
}

__END__
import blastFile='blastx.txt' fastaFile='query.fna' meganFile='blastx.rma' maxMatches=25 minScore=50.0 maxExpected=1.0 topPercent=100.0 minSupport=1 minComplexity=0.30 useSeed=false useCOG=false useKegg=false paired=false useIdentityFilter=false textStoragePolicy=1 blastFormat=BlastX mapping='Taxonomy:BUILT_IN=true';
unCollapse nodes=all;
select nodes=all;
export what=alignment file='example/alignment-%c-%r-%n.fasta' data=Taxonomy classId=selected asConsensus=false useEachReadOnlyOnce=false useEachReferenceOnlyOnce=false refSeqOnly=false translateCDNA=false saveDiversityExtrapolation=false;
quit;
