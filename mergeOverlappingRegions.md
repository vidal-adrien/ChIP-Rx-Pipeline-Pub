# mergeOverlappingRegions.sh

Searches for overlaps of the genomic regions of the between two files and generates merged regions bed file out of selected regions.

Usage: `mergeOverlappingRegions -a <bed/gff/vcf/bam> -b <bed/gff/vcf/bam> [-l <min overlap in bp>]`

Options:
*  `-a`  file A. Any format accepted by bedtools. Each feature in A is compared to B in search of overlaps.
*  `-b` files B. Any format accepted by bedtools. Use -b mutiple times to search overlaps of a in more than one file.
*  `-l` overlap. Value between 0 and 1. Minimum overlap between two regions from A and B respectively as a fraction of either region.

![mergeOverlappingRegions diagram](./images/mergeoverlapdiag.png)

https://github.com/vidal-adrien/PNDS-pipelines/blob/37f2902eeba18af5d5d301a7f2c2201e9ce43b89/bash/mergeOverlappingRegions.sh#L1-L66
