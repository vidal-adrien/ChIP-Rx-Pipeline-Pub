# PNDS Genomics pipelines

## Pipeline guides

Step-by-step guides to realize analyses.

[ChIP-Seq pipeline](chipseq.md)

[ChIP-Rx pipeline](chiprx.md)

## Custom tools 
Scripts used in the above guides.

[bedFromFasta.pl:](bedFromFasta.md): Perl script. Creates a `.bed` table of the full length of the sequences from a `.fasta` file.

[bedFromGff.pl](bedFromGff.md): Perl Script. Creates a `.bed` table of the regions from a `.gff` file. With the possibility to specify which tag contains the ID, to filter by feature type and to enforce ID uniqueness.

[mergeOverlappingRegions.sh](mergeOverlappingRegions.md): Bash script. Searches for overlaps of the genomic regions of the between two files and generates merged regions bed file out of selected regions.

## Genomics resources:

[Col-CEN_blacklist.bed](resources/Col-CEN_blacklist.bed): A blacklist of aberrant regions for the Col-CEN genome.