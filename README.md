pAss - protein-guided Assembler
====

Method for counting genes from a metagenomic survey

## Dependencies

* [MEGAN](http://ab.inf.uni-tuebingen.de/software/megan/)
* [Muscle](https://github.com/Homebrew/homebrew-science). Use `brew` for installation.
* Perl
    * >=v5.20.0 is required:
       * Use [tokuhirom/plenv](https://github.com/tokuhirom/plenv) to install a local version of the latest perl (>=5.21.0)
       * Install **CPANMINUS**, using `plenv install-cpanm`

    * Install dependencies
        * Install Carton using `$ cpanm Carton`.
        * run `carton` in package’s root directory
* R for plots

## Pipeline

The pAss pipeline requires one to provide contigs grouped by the KEGG’s Ortholog groups built using NEWBLER (`runAssembly` >= 2.6 (20110517_1502)).
The scripts should be run in series from 00 to XX.

| Script Name | Description                                                                                                                                    |
| ---          | ---                                                                                                                                            |
| pAss.00     | Given contigs assembled using reads binned into relavant KEGG orthologs we blast them against known reference sequences in the same categories |
| pAss.01     | Calls MEGAN to output the blastx alignment of contigs same KO Refseq sequences                                                                 |
| pAss.03     | Calls the PASS::Alignment package; Details below                                                                                               |
| pAss.04     | Generates diagnostic plots                                                                                                                     |
| pAss.11     | Diagnostic                                                                                                                                     |
| pAss.12     | Outputs as fasta MAX DIVERSITY region for each contig                                                                                          |

## pAss.03 PASS::Alignment

### Input

MSA of contigs to protein reference sequences.

Example 
```
>ref|YP_001105148.1| zinc-containing dehydrogenase [Saccharopolyspora erythraea NRRL
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------W--G--G--L--F--S--D--L--V--R--V--P--W--A--E--A--M--L--V--P--L--P--T--G--P--D--P--V--A--M--A--S--A--S--D--N--------------------------------------------------G--A--R--V--L--V--V--A--R--G--S--I--G--L--Y--V--C--D--I--A--R--A--L--G--A--G--D--V--L--Y--V--D--P--D--P--A--H--R--A--L--A--E--Q--Y--G--A--R--T--A--E--E--I--E--P--
>contig00760  length=176   numreads=21 (rev)
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GGAGCCCGCGTCCTAGTGATGAGCGGCGGCAGCATCGGCCTATATGTCTGTGACATCGCGAGGGCGCTTGGAGCGGCCGAGGTTCTCTACGTCGACCGCGATTCCAGGCGCCGCTCAATTGCGGCCGGCTACGGGGCCAAGACCGCTGAAGCGATTGAGcCA
>contig01620  length=117   numreads=6
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------tggggcggagCTCTCTCGGAGCGCGTGAAAGTTCCGTGGGCCGAAGCGATGCTCCGGCCGATTCCCGCAGGCTTAGATGCCCTGCATTTGTCGAGCCTGAGTGACAAC------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------]
```

### Procedure

1. MUSCLE to generate protein based MSA of reference sequences

2. Aligns and sets coordinates of contig<->REFSEQ(Protein) to the protein MSA from above.

The MAX diversity region is not calculated for a KO if:

1. GAPS (ie. more gaps than are bases)
  * When the the gap ratio (NT:GAP) exceeds 1:10

2. KOs where there are too few contigs are not considered 
  * this or may not be IDEAL cause we’re looking for cases where there’s a large amt of RNA reads mapped to a low number o contigs
