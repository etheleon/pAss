#!/usr/bin/env perl

use FindBin qw/$Bin/;
use Modern::Perl '2015';
use experimental qw/signatures/;
use Parallel::ForkManager;
use Pod::Usage;
use autodie;
use Getopt::Lucid qw/:all/;

my @specs = (
    Param("kodb|d")->default("/export2/home/uesu/db/konr"),
    Param("contigs|c"),
    Param("output|o"),
    Switch("format|f")->default(0),
    Switch ("help|h"),
    Param("threads|t")->default(1)
);

my $opt = Getopt::Lucid->getopt( \@specs );
pod2usage(-verbose=>2) if $opt->get_help;
$opt->validate({'requires' => ['contigs', 'kodb', 'output']});

$|++;
my $threads = $opt->get_threads;
my $pm = Parallel::ForkManager->new($threads);

my $outputdir = $opt->get_output;
system "mkdir $outputdir" unless -d $outputdir;
my $kodb         = $opt->get_kodb;
my %kohash;
my $toFormat = $opt->get_format;


#Step0: Only makeblastdb for KOs which have assemblies completed, ie. with folders in the contig dir
my $contig = $opt->get_contigs;
$kohash{"ko:$_"}++ for map { (split/\//)[-1] } <$contig/*>;

#Step1::Index KOs and build the blast library
#{
    #m/ko:K\d{5}(?!\.phr|\.pin|\.pnd|\.pni|\.pog|\.psd|\.psi|\.psq)/;
    #if(m/ko:K\d{5}$/){
        #my $isKOI = exists $kohash{$&};
        if($toFormat)
        {
            for my $key (keys %kohash) (<"$kodb/*">)
            {
                `makeblastdb -dbtype prot -in $kodb/$& -parse_seqids -out $kodb/$&`;
                say STDERR $key;
            }
        }
    #}
}
say STDERR "    #stored KOs";

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
    `blastx \\
        -db $kodb/$& \\
        -query $contigPath/$1/454AllContigs.fna \\
        -out $output/$1.blastx \\
        -evalue 1e-5 \\
        -seg no \\
        -num_descriptions 10 \\
        -num_alignments 10 \\
        -num_threads $threads`;
    say STDERR "$ko blasted";
}

say STDERR "Blasting is complete";
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


