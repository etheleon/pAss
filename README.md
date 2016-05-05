pAss - protein-guided Assembler
====

[![DOI](https://zenodo.org/badge/19045/etheleon/pAss.svg)](https://zenodo.org/badge/latestdoi/19045/etheleon/pAss)

pAss allows for the enumeration of assembled genes/transcripts within a functional bin.


## Dependencies
* Blast (supports legancy blastall and formatdb)
* [MEGAN5](http://ab.inf.uni-tuebingen.de/software/megan/) + [xvfb-run](http://manpages.ubuntu.com/manpages/lucid/man1/xvfb-run.1.html)
* [Muscle](https://github.com/Homebrew/homebrew-science). Use `brew` for installation.
* Perl
    * v5.2X.X is required:
       * Use [tokuhirom/plenv](https://github.com/tokuhirom/plenv) to install a local version of the latest perl (>=5.21.0)
       * Install **CPANMINUS**, using `plenv install-cpanm`
       * run `./INSTALL`
* R
  * ggplot2
  * dplyr
  * magrittr
  * Biostrings (from Bioconductor)
  * MetamapsDB from github 
* Blast++


## Usage

`$` denotes the terminal prompt

```
$ ./maxDiversity --megan </path/to/MEGAN> --meganLicense <path/to/MEGAN5-academic-license.txt> -f --outputDIR </path/to/output/dir> --contigs ./example/data/contigs/ --refseqKO ./example/refSeqProtDB/ -t 20
```

### Required Input

1. Contigs are assembled from functionally binned reads (eg. NEWBLER 2.6 (20110517_1502))
2. Reference sequences grouped by their gene families
## Description of pipeline

The pAss pipeline requires one to provide contigs grouped by the their Ortholog groups.
The scripts should be run in series from 00 to XX.

| Script Name | Description                                                                                                                                    |
| ---          | ---                                                                                                                                            |
| pAss.00     | Given contigs assembled using reads binned into relavant KEGG orthologs we blast them against known reference sequences in the same categories |
| pAss.01     | Calls MEGAN to output the blastx alignment of contigs same KO Refseq sequences                                                                 |
| pAss.03     | Calls the PASS::Alignment package; Details below                                                                                               |
| pAss.04     | Generates diagnostic plots                                                                                                                     |
| pAss.10     | Scans for maxdiversity region                                                                                                                                      |
| pAss.11     | Outputs MAX diversity sequence                                                                                                                                     |
| pAss.12     | Outputs as fasta MAX DIVERSITY region for each contig                                                                                          |

### pAss.03 PASS::Alignment
This repo is itself a perl module

#### Input

Takes a MSA of contigs to protein reference sequences.

Example
```
>ref|YP_001105148.1| zinc-containing dehydrogenase [Saccharopolyspora erythraea NRRL
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------W--G--G--L--F--S--D--L--V--R--V--P--W--A--E--A--M--L--V--P--L--P--T--G--P--D--P--V--A--M--A--S--A--S--D--N--------------------------------------------------G--A--R--V--L--V--V--A--R--G--S--I--G--L--Y--V--C--D--I--A--R--A--L--G--A--G--D--V--L--Y--V--D--P--D--P--A--H--R--A--L--A--E--Q--Y--G--A--R--T--A--E--E--I--E--P--
>contig00760  length=176   numreads=21 (rev)
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GGAGCCCGCGTCCTAGTGATGAGCGGCGGCAGCATCGGCCTATATGTCTGTGACATCGCGAGGGCGCTTGGAGCGGCCGAGGTTCTCTACGTCGACCGCGATTCCAGGCGCCGCTCAATTGCGGCCGGCTACGGGGCCAAGACCGCTGAAGCGATTGAGcCA
>contig01620  length=117   numreads=6
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------tggggcggagCTCTCTCGGAGCGCGTGAAAGTTCCGTGGGCCGAAGCGATGCTCCGGCCGATTCCCGCAGGCTTAGATGCCCTGCATTTGTCGAGCCTGAGTGACAAC------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------]
```

#### Procedure

1. MUSCLE to generate protein based MSA of reference prot sequences.
2. First Aligns NT contigs to best matched prot REFSEQ then to global REFSEQ protein MSA from above.
3. The MAX diversity region is not calculated for a KO if:

* GAPS (ie. more gaps than are bases)
  * When the the gap ratio (NT:GAP) exceeds 1:10

* KOs where there are too few contigs are not considered 
  * this or may not be IDEAL cause we’re looking for cases where there’s a large amt of RNA reads mapped to a low number o contigs
