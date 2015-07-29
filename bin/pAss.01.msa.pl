#!/usr/bin/env perl

use FindBin qw/$Bin/;
use local::lib "$Bin/../local";
use Modern::Perl '2015';
use autodie;
use Pod::Usage;
use experimental qw/signatures/;
use Getopt::Lucid qw( :all );

my @specs = (
    Param("queryFasta|q"),
    Param("blastFile|b"),
    Param("output|o"),
    Param("megan|m")->default("/export2/home/uesu/local/megan/MEGAN"),
    Param("meganLicense|l")->default("/export2/home/uesu/downloads/MEGAN5-academic-license.txt"),
    Switch ("help|h")
);

my $opt = Getopt::Lucid->getopt( \@specs );
pod2usage(-verbose=>2) if $opt->get_help;
$opt->validate({'requires' => ['queryFasta', 'blastFile', 'output']});

my $license       = $opt->get_meganLicense;
my $query         = $opt->get_queryFasta;
my $blast         = $opt->get_blastFile;
my $out           = $opt->get_output;
$out =~ s/\/$//g;
my $megan         = $opt->get_megan;

die "$query does not exist"    unless -f $query;
die "$blast does not exists"   unless -f $blast;
die "$megan does not exists\n" unless -f $megan;

my $temp = rand().time();
&prep();
&runMegan();

system "rm $temp";

sub runMegan ($count = 1)
{
    if($count <= 20){
        $count++;
    }else{
        say STDERR "$query has failed";
        system "rm $temp";
        exit 1;
    };
    unlink "$out.lock" if -e "$out.lock";
    unlink "$out.log" if -e "$out.log";

    my $scr = int(30000 * rand());

    #try new screen number if its taken;
    $scr = int(30000 * rand()) while -e "/tmp/.X$scr-lock";

    my $signal = `xvfb-run -n $scr -f $out.lock -e $out.log $megan -g -d -E -L $license -c $temp`;
    if($signal =~ m/Writing \d+ reads to file/sm)
    {
        unlink "$out.lock", $temp, "$temp.rma";
        exit 1;
    }else
    {
        say STDERR "reattempt";
        say STDERR $signal;
        runMegan($count);
    }
}

sub prep {
    my $CMD = join '', <DATA>;
    #problem with script
    $CMD =~ s/\|\|blastx\.txt\|\|/$blast/;   #text file
    $CMD =~ s/\|\|blastx\.rma\|\|/$temp.rma/;
    $CMD =~ s/\|\|query\.fna\|\|/$query/;
    $CMD =~ s/\|\|example\|\|/$out/;

    open my $cmdIO, ">", $temp;
    print $cmdIO $CMD;
    close $cmdIO;

    system"mkdir -p $out" unless -d $out;
}

=pod

=head1 NAME

    refMSA - constructing refseq MSA

=head1 SYPNOSIS

    Aligns contig sequences with ortholog group reference sequences using megan

=head1 OPTIONS

=over 4

=item --queryFasta -q

    the fastaFile containing contigs

=item --blastFile -b

    the output from blasting contigs against reference sequences

=item --megan -m

    path to megan executable
=item --meganLicense -l

    path to MEGAN's academic license

=back
=cut

__END__
import blastFile='||blastx.txt||' fastaFile='||query.fna||' meganFile='||blastx.rma||' maxMatches=25 minScore=50.0 maxExpected=1.0 topPercent=100.0 minSupport=1 minComplexity=0.30 useSeed=false useCOG=false useKegg=false paired=false useIdentityFilter=false textStoragePolicy=1 blastFormat=BlastX mapping='Taxonomy:BUILT_IN=true';
unCollapse nodes=all;
select nodes=all;
export what=alignment file='||example||/alignment-%c-%r-%n.fasta' data=Taxonomy classId=selected asConsensus=false useEachReadOnlyOnce=false useEachReferenceOnlyOnce=false refSeqOnly=false translateCDNA=false saveDiversityExtrapolation=false;
quit;
