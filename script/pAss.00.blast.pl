#!/usr/bin/env perl

use Modern::Perl '2015';
use experimental qw/signatures/;
use Parallel::ForkManager;
use Pod::Usage;
use autodie;
use Getopt::Lucid qw(:all);

my @specs = (
    Param("kodb|d")->default("/export2/home/uesu/db/konr"),
#    Param("threads|t")->default(1),
#    Param("batch|b")->default("pAss00_BATCH"),
    Param("contigs|c")->default("/export2/home/uesu/reDiamond/out/assm.0200.newbler"),
    Param("output|o")->default("/export2/home/uesu/reDiamond/out/assm.0300"),
    Switch("format|f")->default(0),
    Switch ("help|h"),
    Param("threads|t")->default(1)
);

my $opt = Getopt::Lucid->getopt( \@specs );
pod2usage(-verbose=>2) if $opt->get_help;

$|++;

my $pm = Parallel::ForkManager->new($opt->get_threads);

my $kodb         = $opt->get_kodb;
my %kohash;

#Step1::Index KOs and build the blast library
for (<"$kodb/*">)
{
    m/ko:K\d{5}/;
    $kohash{$&}++;
    if($opt->get_format){
        system "formatdb -i $kodb/$& -o T - $kodb/$&";
        say STDERR $_;
    }
}
say STDERR "stored KOs";

#open my $batchOutput, ">", $opt->get_batch;
for my $indivKO (keys %kohash)
{

    $pm->start and next;
    runBLAST($indivKO, $opt->get_contigs, $opt->get_output);
    $pm->finish;
}

$pm->wait_all_children;


sub runBLAST($ko, $contigPath, $output)
{
    $ko =~ m/ko:(K\d{5})/;
    unless (-e "$contigPath/$1/454AllContigs.fna")
    {
        return;
    }
    system "blastall -p blastx -v 10 -b 10 -F F -e 1e-5 -a 4 -d $kodb/$& -i $contigPath/$1/454AllContigs.fna -o $output/$1.blastx"
}

=pod

=head1 NAME

    refBLAST - Blasting KEGG Orthology binned contigs against reference genomes

=head1 SYPNOSIS

    Runs blast of contigs binned int he same functional group, KOs, against reference sequences annoted from the same KEGG Orthology Groups
    $ refBLAST -d /export2/home/uesu/db/konr -c out/assm.0200.newbler -o out/pAss00 -f

=head1 OPTIONS

=over 4

=item --kodb -d

    path to the root directory housing reference sequences for each KO in separate children directories

=item --batch -b

    path to directory to store the commonds to be executed by GNU parallel or submitted as a job in a SGE system

=item --contigs -c

    path to the contigs, parent dir storiing child directories each from the output form running the assembly process (continue from the pAss pipeline)

=item --output -o

    path to store the output of the blast

=item --format -f

    to run formatdb first

=back

=cut


