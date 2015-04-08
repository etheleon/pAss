pAss - protein-guided Assembler
====

Method for counting genes from a metagenomic survey

## Dependencies
* MEGAN
* Muscle 
* Perl 
 * Set::IntervalTree (GCC > 4.1.2)
 * Moo
 * namespace::clean
 * Bio::SeqIO
 * Statistics::Basic
 * List::MoreUtils
 * Data::Dumper

## Alignment
The MAX diversity region is not calculated for a KO if:

1. GAPS (ie. more gaps than are bases)
  * When the the gap ratio (NT:GAP) exceeds 1:10
2. KOs where there are too few contigs are not considered 
  * this or may not be IDEAL cause we’re looking for cases where there’s a large amt of RNA reads mapped to a low number o contigs
