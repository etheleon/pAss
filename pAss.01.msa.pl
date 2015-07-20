#!/usr/bin/env perl

use Modern::Perl '2015';
use autodie;
use Pod::Usage;
use experimental qw/signatures/;
use Getopt::Lucid qw( :all );

my @specs =(
    Param("queryFasta|q"),
    Param("blastFile|b"),
    Param("output|o"),
    Param("megan|m"),
    Switch ("help|h")
);

my $opt = Getopt::Lucid->getopt( \@specs );
pod2usage(-verbose=>2) if $opt->get_help;
$opt->validate({'requires' => ['queryFasta', 'blastFile', 'output', 'megan']});

my $query = $opt->get_queryFasta;
my $blast = $opt->get_blastFile;
my $out = $opt->get_output;
my $megan = $opt->get_megan;

die unless -f $query;
die unless -f $blast;
die unless -f $megan;

my $temp = rand().time();

&prep();
&runMegan();

sub runMegan ($count = 1)
{
    if($count <= 20){
        $count++;
    }else{
        say STDERR "$query has failed";
        exit 1;
    };
    unlink "$out.lock" if -e "$out.lock";
    unlink "$out.log" if -e "$out.log";

    my $scr = int(30000 * rand());

    #try new screen number if its taken;
    $scr = int(30000 * rand()) while -e "/tmp/.X$scr-lock";

    my $signal = `xvfb-run -n $scr -f $out.lock -e $out.log $megan -g -d -E  -c $temp`;
    #cant check for xvfb-run's own error [xc's version]
    if($signal =~ m/Writing \d+ reads to file/sm)
    {
        unlink "$out.lock", $temp, "$temp.rma";
        exit 1;
    }else{
        say STDERR "reattempt";
        say STDERR $signal;
        runMegan($count);
    }
}

sub prep {
    my $CMD = join '', <DATA>;
    $CMD =~ s/blastx.txt/$blast/;   #text file
    $CMD =~ s/blastx.rma/$temp.rma/;
    $CMD =~ s/query.fna/$query/;
    $CMD =~ s/example/$out/;

    open my $cmdIO, ">", $temp;
    print $cmdIO $CMD;
    close $cmdIO;

    mkdir $out unless -d $out;
}

=pod

=head1 NAME

    refMSA

=head1 Aligns contig sequences with ortholog group reference sequences using megan

=head1 OPTIONS

=over 4

=item --queryFasta -q

    the fastaFile containing contigs

=item --blastFile -b

    the output from blasting contigs against reference sequences

=item --megan -m

    path to megan executable

=back
=cut

__END__
import blastFile='blastx.txt' fastaFile='query.fna' meganFile='blastx.rma' maxMatches=25 minScore=50.0 maxExpected=1.0 topPercent=100.0 minSupport=1 minComplexity=0.30 useSeed=false useCOG=false useKegg=false paired=false useIdentityFilter=false textStoragePolicy=1 blastFormat=BlastX mapping='Taxonomy:BUILT_IN=true';
unCollapse nodes=all;
select nodes=all;
export what=alignment file='example/alignment-%c-%r-%n.fasta' data=Taxonomy classId=selected asConsensus=false useEachReadOnlyOnce=false useEachReferenceOnlyOnce=false refSeqOnly=false translateCDNA=false saveDiversityExtrapolation=false;
quit;
