use strict;
use warnings;

use Modern::Perl '2015';
use experimental qw/postderef signatures/;
use PASS::Alignment;
use Test::More;
use FindBin;
use Data::Dumper;
my $module = 'PASS::Alignment';

use_ok($module)
    or BAIL_OUT("Does $module.pm compile?  Does it end with 1; ?");
can_ok($module, 'storeRefseq')
    or BAIL_OUT("Missing package $module; or missing sub storeRefseq()");
can_ok($module, 'assignContig2ref')
    or BAIL_OUT("Missing package $module; or missing sub assignContig2ref()");
can_ok($module, 'runMuscle')
    or BAIL_OUT("Missing package $module; or missing sub runMuscle()");
can_ok($module, 'buildContigMSA')
    or BAIL_OUT("Missing package $module; or missing sub buildContigMSA()");
can_ok($module, 'grepRefSeqID')
    or BAIL_OUT("Missing package $module; or missing sub grepRefSeqID()");
can_ok($module, 'minmax')
    or BAIL_OUT("Missing package $module; or missing sub minmax()");

my $alignment = PASS::Alignment->new(
    refseqFasta     =>  "$FindBin::Bin/data/refSeqProtDB/ko\:K00001", #reference sequences
    alignmentFiles  =>  [glob("$FindBin::Bin/data/pAss01/K00001/alignment-*")],
    outputPrefix    =>  "out",
)                                                                              ;

$alignment->storeRefseq();

ok scalar keys $alignment->{refseq}->%* == 772, "Loading: read allContigs In";
my $sub = $module->can('hey');

$alignment->assignContig2ref();
my @missing;
for my $contigKey (keys $alignment->{contigs}->%*)
{
    for my $ref ($alignment->{contigs}->{$contigKey}->{parentREF})
    {
        if ($ref eq "")
        {
            push @missing, $contigKey;
        }
    }
}

is scalar @missing, 0, "Everyone is partnered";



done_testing();

__DATA__
