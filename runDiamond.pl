#!/usr/bin/env perl


use strict;
use warnings;

die "$0 <#threads> <diamond DB without fileExtention> <fastQFILE> <scratch> <scratchOutputFile> <nasOutputDir>
    eg. runDiamond.pl 22 /scratch/uesu/db/nrfull /export2/ulu_pandan2011/data/batch1_gDNA_illumina/filtered_fastq/s_1_1.filtered.fastq /scratch/uesu/ /scratch/uesu/db/out/s_1_1.filtered.m8 /export2/home/uesu/diamonoutput\n" unless $#ARGV == 5;


my ($threads, $nrDB, $fq, $scratch, $outputFile, $storageDir) = @ARGV;
mkdir $storageDir unless -d $storageDir;
#print "diamond blastx -v -p $threads -d $nrDB -q $fq -o $outputFile -t $scratch --index-chunks 1", "\n";
system "/export2/home/uesu/local/bin/diamond blastx -v -p $threads -d $nrDB -q $fq -o $outputFile -t $scratch --index-chunks 2";
my $filename = (split /\//, $outputFile)[-1];
#print "mv $outputFile $storageDir/$filename", "\n";
system "mv $outputFile $storageDir/$filename";

__END__
qsub -N jobs_1_1.filtered -m abe -M wesley@bic.nus.edu.sg -b y -V -cwd -pe smp 22 script/runDiamond.pl <ur arguments>
