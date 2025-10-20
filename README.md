# Plant Nuclear Dynamics & Signaling (PNDS) pipelines
[![DOI](https://zenodo.org/badge/633768138.svg)](https://doi.org/10.5281/zenodo.17397579)

This is a public collection of documentation about ChIP-seq and ChIP-Rx pipelines used by the PNDS (Plant Nuclear Dynamics & Signaling) Team led by Clara Bourbousse & Fredy Barneche.

Scripts and documentation by Adrien Vidal.

## Pipeline guides

Step-by-step guides to genomic analysis pipelines used by the team.

[ChIP-Seq pipeline](chipseq.md)

[ChIP-Rx pipeline](chiprx.md)

## Custom tools
Scripts used in the above guides.

[bedFromFasta.pl](bedFromFasta.md): Perl script. Creates a `.bed` table of the full length of the sequences from a `.fasta` file.

[bedFromGff.pl](bedFromGff.md): Perl Script. Creates a `.bed` table of the regions from a `.gff` file. With the possibility to specify which tag contains the ID, to filter by feature type and to enforce ID uniqueness.

[mergeOverlappingRegions.sh](mergeOverlappingRegions.md): Bash script. Uses a comination of bedtools functions to search for overlaps of the genomic regions of the between two files and generate merged regions bed file out of selected regions.

## Genomics resources

Genomic resources used by the the team when applying the above pipelines to *Arabidopsis thaliana* experiments. 

**Reference genomes:**
*  [⇗TAIR10_chr_all.fas.gz](https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz): *Arabidopsis thaliana* TAIR10 genome assembly.
*  [⇗Col-CEN_v1.2.fasta.gz](https://github.com/schatzlab/Col-CEN/blob/main/v1.2/Col-CEN_v1.2.fasta.gz): *Arabidopsis thaliana* Col-CEN genome assembly.
*  [⇗Download page](https://www.ebi.ac.uk/ena/browser/view/GCA_028009825) for *Arabidopsis thaliana* Col-CC genome assembly.
*  [⇗Download page](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001215.4/) for *Drosophila melanogaster* release 6 genome assembly.

**Annotation:**
*  [Araport11_GFF3.gene.201606.bed](resources/Araport11_GFF3.gene.201606.bed): Araport 11 annotation for genes on *Arabidopsis thaliana* TAIR10 genome assembly as a `.bed` file. Converted from the [june 2016 annotation](https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/archived/Araport11_GFF3_genes_transposons.Jun2016.gff.gz).
*  [Araport11_GFF3.TAIR10.transposable_element.201606.bed](resources/Araport11_GFF3.TAIR10.transposable_element.201606.bed): Araport 11 annotation for transposable elements on *Arabidopsis thaliana* TAIR10 genome assembly. Converted from the [june 2016 annotation](https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/archived/Araport11_GFF3_genes_transposons.Jun2016.gff.gz).
*  [Col-CEN_v1.2_genes.araport11.gene.bed](resources/Col-CEN_v1.2_genes.araport11.gene.bed): Lifted Araport 11 annotation for genes on *Arabidopsis thaliana* Col-CEN genome assembly. Converted from the [`.gff3` annotation](https://www.arabidopsis.org/download_files/Genes/Col-CEN%20genome%20assembly%20release/ColCEN_GENES_Araport11.gff3.gz).

**Blacklists:**
*  [TAIR10_blacklist.bed](resources/TAIR10_blacklist.bed): A blacklist of aberrant regions for the *Arabidopsis thaliana* TAIR10 genome.
*  [Col-CEN_blacklist.bed](resources/Col-CEN_blacklist.bed): A blacklist of aberrant regions for the *Arabidopsis thaliana* Col-CEN genome.
