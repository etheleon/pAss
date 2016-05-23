PADI - Protein-guided Assembly and Diversity Indexing
====

[![DOI](https://zenodo.org/badge/19045/etheleon/pAss.svg)](https://zenodo.org/badge/latestdoi/19045/etheleon/pAss)

## Description

PADI enumerates and analyses protein guided assemblies of proteinaceous functional bins. 
We look to the regions of sequence maximal diversity determined by global multiple sequence alignment (MSA)

### Procedure

1. BlastX of nucleotide contigs against reference prot sequences
2. Assign NTcontigs to best matched prot REFSEQ2.
3. MUSCLE generates MSA of reference prot sequences
4. Align NTcontigs to global REFSEQ MSA
3. Search for MAX diversity region

#### Filters

* GAPS (ie. more gaps than are bases)
  * When the the gap ratio (NT:GAP) exceeds 1:10

* KOs where there are too few contigs are not considered
  * this or may not be IDEAL cause we’re looking for cases where there’s a large amt of RNA reads mapped to a low number o contigs


## Installation

### Docker Installation (Recommended)

[Install](https://docs.docker.com/engine/installation/) Docker

#### Usage: SINGLE COPY GENES

```
$ cp SingleCopyGene SCG
$ mkdir SCG/out SCG/misc
$ #Place MEGAN5 license file inside misc
$ cp MEGAN5-academic-license.txt SCG/misc/
```

```
docker run --rm \
    -v `pwd`/SCG/data/konr:/data/refSeqProtDB \
    -v `pwd`/SCG/data/newbler:/data/contigs \
    -v `pwd`/SCG/out:/data/out \
    -v `pwd`/SCG/misc:/data/misc \
    pass \
    /tmp/pAss/maxDiversity --outputDIR /data/out --format --threads 8 --refseqKO /data/refSeqProtDB  --contigs /data/contigs  --megan /usr/local/bin/MEGAN --meganLicense /data/misc/MEGAN5-academic-license.txt 
```

### Source Installation 

#### Dependencies/Pre-requisites

* Tools
    * [MEGAN5](http://ab.inf.uni-tuebingen.de/software/megan/) 
    * [xvfb-run](http://manpages.ubuntu.com/manpages/lucid/man1/xvfb-run.1.html)
    * [Muscle](https://github.com/Homebrew/homebrew-science). Use `brew` for installation.
    * Blast+

* Perl
    * v5.2X.X is required:
       * Use [tokuhirom/plenv](https://github.com/tokuhirom/plenv) to install a local version of the latest perl (>=5.21.0)
       * Install **CPANMINUS**, using `plenv install-cpanm`

* R
  * ggplot2
  * dplyr
  * Biostrings (from Bioconductor)


#### Usage

```
$ ./maxDiversity --megan </path/to/MEGAN> --meganLicense <path/to/MEGAN5-academic-license.txt> -f --outputDIR </path/to/output/dir> --contigs ./example/data/contigs/ --refseqKO ./example/refSeqProtDB/ --theads 20
```

### Required Input


1. Contigs are assembled from functionally binned reads (eg. NEWBLER 2.6 (20110517_1502))
2. Reference sequences grouped by their gene families

Use [pipeline](https://github.com/quanyu2015/ngs_pipeline) to get from raw reads to end of procedure.

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

## Publication

In preparation

## Future

Include a sister software which uses protein HMMs on top of sequence similarity as a method to search for distantly related sequences.


