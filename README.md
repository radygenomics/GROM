# GROM
Version 1.0.1




## Description, Installation, and Usage


### DESCRIPTION

GROM detects a comprehensive range of variants, SNVs, indels (deletion and insertion), structural variants (deletion, duplication, insertion, inversion, translocation), and CNVs (deletion and duplication).  GROM integrates paired-end, split-read, and read depth information.  GROM is optimally designed for whole genome sequencing (Illumina paired-end or single-end) data aligned with BWA mem and has been used to analyze exome and RNA-seq data.  GROM analysis is limited to SNVs and indels for exome and RNA-seq data.  GROM has been tested with BWA, however we recommend using BWA mem since GROM utilizes BWA mem's split-read mapping information.  When analyzing single-end data, GROM's insertion, inversion, and translocation output will be very limited.  This is a common limitation of single-end sequencing analysis not limited to GROM.  GROM is capable of duplicate read filtering.  Duplicate read filtering is similar to Picard's MarkDuplicates with an exception that GROM does not consider soft-clipping of a read's mate.  This is due to the speed optimization of GROM ("thunder" in Russian), and GROM is "lightning" fast completing analysis many times faster than leading algorithms. When run single-threaded, GROM requires approximately 13 GB RAM for a human genome.  Running GROM with 24-threads requires 128 GB RAM for a human genome.

GROM outputs variant predictions in VCF-format.  Translocations are output to a separate file with the extension ".ctx.vcf" instead of ".vcf".

GROM utilizes SAMtools (version 1.3.1 included with the GROM distribution) for accessing BAM files.  GROM has been successfully compiled and tested with:

- CentOS 6.8; gcc 4.4.7
- Ubuntu 14.04; gcc 4.4.8
- Ubuntu 14.04; gcc 6.3





### INSTALLATION

The easiest way to use GROM is to use the precompiled binary (created using gcc 4.4.7 on CentOS 6.8) in ./dist/

If desired, GROM may be created from source code by running the following from the command line:

```
make
```





### USAGE

GROM is intended to be easy to use. For typical usage, such as with a human genome, GROM requires 3 input parameters:

1) BAM file
2) Reference file (in fasta format)
3) Output file

Note, the BAM file must be sorted by coordinate and indexed.  

Example using the test data with the precompiled GROM:

```
./dist/GROM -i ./test_data/tilapia_SAMD00023995_GL831235-1.bam -r ./test_data/oreNil2_GL831235-1.fa -o tilapia_SAMD00023995_GL831235-1_predictions.vcf
```

or if compiling from source code:

```
./bin/GROM -i ./test_data/tilapia_SAMD00023995_GL831235-1.bam -r ./test_data/oreNil2_GL831235-1.fa -o tilapia_SAMD00023995_GL831235-1_predictions.vcf
```

Output file should be similar to results in ./test_data/test_outuput_tilapia_SAMD00023995_GL831235-1.vcf and ./test_data/test_outuput_tilapia_SAMD00023995_GL831235-1.ctx.vcf


We recommend using the '-M' parameter to filter duplicate reads.  This will likely improve the variant prediction ability of GROM and often result in faster run times. Example:

```
./dist/GROM -M -i ./test_data/tilapia_SAMD00023995_GL831235-1.bam -r ./test_data/oreNil2_GL831235-1.fa -o tilapia_SAMD00023995_GL831235-1_predictions.vcf
```


Default is for a female genome.  For a male genome, use '-g 1'.  Example:

```
./dist/GROM -g 1 -i [BAM file] -r [Reference file] -o [Output file]
```


For genomes with ploidy not equal to 2, use '-p [ploidy].  For instance, if the ploidy is 4:

```
./dist/GROM -p 4 -i [BAM file] -r [Reference file] -o [Output file]
```


For multi-threading, use -P [processes].  [processes] indicates number of processes (chromosomes) to run simultaneously.  Each process requires 2 threads.  For instance, "-P 12" will run 12 processes (12 chromosomes) in parallel and require 24 threads.  Default is single-threaded.  Example of 24-threaded GROM:

```
./dist/GROM -P 12 -i [BAM file] -r [Reference file] -o [Output file]
```


### Other parameters
The remainder of this file describes other command line parameters.  For most cases, the following parameters do not need to be adjusted.


`-b base phred quality score threshold [default=20]`

Bases with phred quality score below the threshold are ignored.


`-q low mapping quality threshold [default=20]`

Reference bases with average read mapping quality below the threshold are analyzed separately.  For CNV detection, GROM creates separate GC distributions for reference bases in low quality regions and high quality regions. For other variants, GROM ignores low quality reads, except for providing read depth information used in probability statistics.


`-v probability score threshold [default=0.001]`

Variants with probability score above the threshold are filtered out.  Not applicable for insertion and CNV detection.  See "-e" for insertion detection and "-V" for CNV detection.


`-e probability score threshold for insertions [default=0.0000000001]`

Insertions with probability score above the threshold are filtered out.


`-V probability score threshold for CNVs [default=0.000000001]`

CNVs with probability score above the threshold are filtered out.


`-A sampling rate [default=2]`

Only affects CNV detection. Controls window sampling rate.  For each window size, a sample distribution is created for each chromosome.  'sampling rate' indicates the number of windows sampled per maximum window size.  For instance if the sampling rate is 2 and the maximum window size is 10,000, then 2 sample windows will be taken for every 10,000 reference bases.

When analyzing small genomes (chromsomes <10,000,000 bases), it may be beneficial to increase the sampling rate to 3 or 4.  CNV detection may be unreliable for chromsomes <1,000,000 bases.


`-B Chromosome maximum length [default=300000000]`

If longest chromosome is longer than 'Chromosome maximum length', increase 'Chromosome maximum length', although this will increase the RAM memory needed to run GROM.


`-D Dinucleotide repeat minimum length [default=20]`

Only affects CNV detection.  Minimum length in bases (not repeats) for dinucleotide repeats.


`-E Dinucleotide repeat minimum standard deviation [default=1.5]`

Only affects CNV detection.  A Dinucleotide repeat must have an average read depth more than 1.5 standard deviations lower than the chromosome average read depth, and vice versa (chromosome average read depth more than 1.5 standard deviations greater than the dinucleotide repeat average read depth).  When conditions are met, a coverage distribution is created for the dinucleotide repeat.  In regions with the dinucleotide repeat, the dinucleotide repeat coverage distribution is used in replace of the weighted GC coverage distribution.


`-K 0=no ranking, 1=use ranks [default=1]`

Only affects CNV detection. If "1", use ranking for variance.


`-L Duplication coverage threshold [default=2]`

Only affects CNV detection.  Read coverage greater than 'Duplication coverage threshold' * (times) average chromosome read depth is adjusted to 'Duplication coverage threshold' * (times) average chromosome read depth.  For single ploidy genomes, this parameter may need to be increased to 3.


`-M`

Turn on duplicate read filtering. [default=OFF]


`-S`

Turns off split-read analysis. [default=ON]


`-U Excessive coverage threshold [default=2]`

Only affects CNV detection.  Block is labelled as excessive coverage when read depth is greater than 'excessive coverage threshold' * chromosome average read depth.  For single ploidy genomes, this paramater may need to be increased to 3.


`-W Window minimum size [default=100]`

Only affects CNV detection.  Smallest window used for detecting CNVs.  For the CNV search, windows sizes are increased until 'Window maximum size' is reached.  If no CNV has been found, the search moves to the next suitable base.  If a CNV was found, the maximum window size is slid to extend the CNV.  No need to adjust this parameter.  


`-X Window maximum size [default=10000]`

Only affects CNV detection.  Largest window used for detecting CNVs.  See description for 'Window minimum size'.


`-Y Minimum number of blocks [default=4]`

Only affects CNV detection.  Number of clustered excessive coverage blocks needed for a genome region to be masked.  GROM does two passes through the genome to detect CNVs.  The first pass masks excessive coverage regions.  The second does not.  A union set of CNVs is created from the two pass throughs.


`-Z Block size [default=10000]`

Only affects CNV detection.  Size of blocks for excessive coverage masking.  If 'Minimum number of blocks' is 4, then a masked region will be at least 40000 (4*10000) bases.


`-a Minimum SNV ratio [default=0.2]`


Only affects SNV detection.  Minimum ratio of SNV to read coverage.


`-d Minimum reads supporting breakpoint [default=2]`

Does not affect SNV detection.


`-h (no argument)`

Display help.


`-j Minimum SV ratio [default=0.05]`

Only affects SV detection.  Minimum ratio of SV evidence to read coverage.


`-k Maximum homopolymer threshold [default=10]`

Indels adjacent to homopolymers exceeding threshold are filtered.


`-m Minimum indel ratio [default=0.125]`

Only affects indel detection.  Minimum ratio of indel evidence to read coverage.


`-n Minimum SNV reads supporting variant [default=3]`

Only affects SNV detection.


`-s Minimum standard deviations for discordance [default=3]`

Determines insert size range for concordant pairs.


`-u Maximum evidence ratio [default=0.25]`

Only affects SV detection.  Weak evidence (soft-clipping, reads with unmapped mates) are included as evidence if ratio of weak evidence to strong evidence (discordant pairs, split-reads) does not exceed "maximum evidence ratio".


`-w Maximum insertion breakpoint threshold [default=10]`

Only affects insertions.  Insertions filtered if start and end breakpoints (on reference) differ by more than threshold.


`-x Minimum average base quality threshold [default=15]`

Only affects SNV detection.  SNVs with average base quality less than threshold are filtered.  Low mapping quality bases are included in average.


`-y Maximum unmapped gap or overlap for split-reads [default=20]`

Split-reads with gap or overlap greater than the threshold are not considered evidence of a variant.


`-z Minimum split-read mapped length (each split) [default=30]`

Split-reads less than the minimum may not support variants based on the split.  Split-reads may still support variants based on discordant pairs (paired-reads).

