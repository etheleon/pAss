use v5.10;
my %hit;
my %msa2ref;
my $ref;
my %pos;
my $seqFull = 'GHSPALGLGTWKSSPQVVGQAVEQALDLGYRHLDCAAIYGNEAEIGATLANAFTKGVVKREELWITSKLWSNAHHPDAVLPALEKTLQDLGLDYLDLYLIHWPVVIQPDVGFPESGDQLLPFTPASLEGTWQALEKAVDLGLCHHIGVSNFSLKKLEMVLSMARIPPAVNQVELHPYLQQSDLLTFANSQNILLTAYSPLGSGDRPAAFQQAAEPKLLTDPVINGIAAEQGCSAAQVLLAWAIQRGTVTIPKSVNPERLEQNLRAADITLTDSEMAKIALLDRHYRYVSGDFWTMPGSPYTLQNLWDE';
print STDERR "this is the file: $ARGV[0]\n";

open(IN, $ARGV[0]) or die;
$/ = "\n>";
while(<IN>)
{
    s/^\s+|\s+$//g;
    s/^\s*>|>\s*$//g;
    s/^\s+|\s+$//g;
#say "\n\nthis is the line: $_";
    #for each sequence
    if(m/^(\S+).*?\n(.+)/s)
    {
        my $id = $1;
        my $seq = $2;
        $seq = uc $seq;
        $seq =~ s/\s+//g;
        my $pos = 0;

        if($id =~ m/^ref\|/)
        ##################################################
        {
            $ref = $id;
            #next unless $good_ref{$ref};

            my $last_end = 0;
#                                                              *
#---------------------------------------------------------------P--A--L--G--L--G--T--W--K--S--S--P--Q--V--V--G--Q--A--V--E--Q--A--L--D--L--G--Y--R--H--L--D--C--A--A--I--Y--G--N--E--A--E--I--G--A--T--L--A--N--A--F--T--K--G--V--V--K--R--E--E--L--W--I--T--S--K--L--W--S--N--A--H--H--P--D--A--V--L--P--A--L--E--K--T--L--Q--D--L--G--L--D--Y--L--D--L--Y--L--I--H--W--P--V--V--I--Q--P--D--V--G--F--P--E--S--G--D--Q--L--L--P--F--T--P--A--S--L--E--G--T--W--Q--A--L--E--K--A--V--D--L--G--L--C--H--H--I--G--V--S--N--F--S--L--K--K--L--E--M--V--L--S--M--A--R--I--P--P--A--V--N--Q--V--E--L--H--P--Y--L--Q--Q--S--D--L--L--T--F--A--N--S--Q--N--I--L--L--T--A--Y--S--P--L--G--S--G--D--R--P--A--A--F--Q--Q--A--A--E--P--K--L--L--T--D--P--V--I--N--G--I--A--A--E--Q--G--C--S--A--A--Q--V--L--L--A--W--A--I--Q--R--G--T--V--T--I--P--K--S--V--N--P--E--R--L--E--Q--N--L--R--A--A--D--I-----T--L--T--D--S--E--M--A--K--I--A--L--L--D--R--H--Y--R--Y--V--S--G--D--F--W--T--M--P--G--S--P--Y--T--L--Q--N--L--W--D--E--
            while($seq =~ m/(([^-]--)+)/g)
            {
                my $start = $-[0];  #the place before the 1st amino acid
                say "FRONT GAPS:",substr($seq, 0, $-[0]);   #until *
                my $fullseq = $1;   #the sequence after the gaps
                say "this is sequence: ", $fullseq;
                (my $stripseq = $fullseq) =~ s/-//g;
                say "stripped till GAP: \n$stripseq\n";
                say "reference sequence: \n", substr($seqFull, $last_end);
                my $index = index(substr($seqFull, $last_end), $stripseq);  #sometimes it isnt a full alignment with the reference sequernce, so i need to know where it begins
                say "\n\nINDEX: $index\n\n";

                my $i = 0;
                while($fullseq =~ m/[^-]/g) #move stepwise from 1 aa to the next aa
                {
                    $i++;
                    say "this is the location: ", $-[0];
                    #say "amino acid $i: ",substr($seq, $-[0], $-[0]+1);
                    my $MSAaaLOC = $start + $-[0] + 1;
                    say "this is MSAaaLOC: ",$MSAaaLOC; #the location of the aa on the msa;
                    $msa2ref{$MSAaaLOC} = $i + $last_end + $index;
                    say "this is the AAloc: ", eval{ $i + $last_end + $index};
                                               #print "\t\t", substr($seq{$id}, $msa2ref{$j}-1, 1), "=", substr($seq, $j-1, 1), "\n";
                }
                $last_end += length $stripseq;
                say "lastEND:",$last_end;
            }
        }
        ##################################################
    }
}

say "$_\t$msa2ref{$_}" for sort {$a <=> $b}keys %msa2ref;
