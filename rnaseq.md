# RNA-Seq pipeline

A pipeline for RNA-Seq analysis.

The following files may be needed for this pipeline:

* The reads sequencing data in `.fastq` format for each sample. May be compressed as a `.gz` archive.
  * For each experimental condition, both an input (sequencing before immuno-precipitation) and IP (sequencing after immuno-precipitation) are needed.
  * In the case of pair-end sequencing, two files per sample with matching read IDs are needed.
* The sequences of potential contaminants to control for in `.fasta` format (optional).
* The sequences of the sequencing adapters in `.fasta` format.
* The sequence of the reference genome in `.fasta` format.
* Optionally, a `.bed` file containing a blacklist of genomic regions to exclude from the map.
* `.bed` files of the genes to analyse. The 4th column of those bed files must contain the identifiers of the regions. Such a file can be produced from a `.gff` annotation using the [bedFromGff.pl](bedFromGff.md) script.

Unless specified otherwise, all code examples are in the bash Unix Shell command language.

The following programs are used in this pipeline:

* [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (optional).
* [fastq_screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/) (optional).
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).
* [STAR](https://github.com/alexdobin/STAR).
* [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (optional).
* The [sambamba](http://lomereiter.github.io/sambamba/) toolkit.
* The [samtools](https://www.htslib.org/) toolkit (optional*).
* The [bedtools](https://bedtools.readthedocs.io/en/latest/index.html) toolkit.
* [MACS2](https://pypi.org/project/MACS2/) (requires python3).
* The [deeptools](https://deeptools.readthedocs.io/en/develop/) toolkit (requires python3).
* The [R](https://www.r-project.org/) language with the following packages:
  * [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
  * [gplots](https://cran.r-project.org/web/packages/gplots/)
  * [ggplot2](https://ggplot2.tidyverse.org/)
  * [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html)

**\*** *All functions performed by samtools can be performed by sambamba. Only one steps require the specific use of sambamba. However, due to having encountered file corruption issues with samtools, sambamba has come to replace it in all tasks that it previously was used for. Samtools may still be used if prefered.*

The following key terms are important to understand this page's instructions: 
  * **Sample**: The set of sequenced reads corresponding to a single experimental sample. For experiments using pair-end sequencing, there will be two sets of reads.
  * **Experimental condition**: Two samples, one taken before immunoprecipitation (input) and the other after (IP) constitute the data for a single experimental condition.
  * **Biological replicate**: Biological replicates are experimental conditions meant to be identical with the same genotype and environmental variables.
  * **Replicate group**: The set of all experimental conditions which are biological replicates of each other.

## <a id="indexing">1) Reference genomes processing and indexing</a>

The referrence genome must be indexed. [STAR](https://github.com/alexdobin/STAR) with the `--runMode genomeGenerate` option will generate an index that it can use in its regular run mode. This step is only needed for the very first analysis made with this reference genome. The index files created by this tool can be used for any subsequent analyses.

The STAR manual advises using the following formula to calculate the value for the `--genomeSAindexNbases` argument of its indexing mode:
$$\min(14,~\frac{1}{2}*log_{2}(\sum_{i=1}^{C}L_i)-1)$$

Where $C$ is the number of chromosomes on the genome assembly to index and $L_i$ is the length in base pairs of a chromosome $i$.

To compute this, a `.bed` file describing the entire genome as genomic regions which can be produced using the [bedFromFasta.pl](bedFromFasta.md) script.

```shell
bedFromFasta.pl -i reference_genome.fasta -o reference_genome.bed
```

This file will also be required for other procedures. 

Here, executing some R code using Rscript allows to more easily compute the formula.

```shell
NBASES=$(Rscript -e "writeLines(as.character(as.integer(
    min(14, log2(sum( 
        read.table('reference_genome.bed')[[3L]]
    ))/2 - 1)
)))")
```

Then the index is built for the STAR aligner using `--runMode genomeGenerate`. The directory given as the `--genomeDir` must be empty.
```shell
STAR --runMode genomeGenerate \
  --runThreadN $THREADS \
  --genomeSAindexNbases NBASES \
  --genomeFastaFiles reference_genome.fasta \
  --genomeDir STAR_index/
```

<a id="bt2indexing"></a>
Finally, for the purpose of performing [contamination screening](#fastqscreen) with [fastq_screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/), the the reference genome should be indexed for [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) using the [bowtie2-build](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer) command.

```shell 
bowtie2-build -threads $THREADS reference_genome.fasta reference_genome
```

The last argument (in this case `"reference_genome"` is the prefix that was used to name each file in the index. This will by default equal the name of the `.fasta` file without the extension if unspecified. Specifying it is useful to differentiate different index builds or specify the path of the directory to contain the index files.

This procedure will produce 6 files respectively to the prefix:

```
reference_genome.1.bt2
reference_genome.2.bt2
reference_genome.3.bt2
reference_genome.4.bt2
reference_genome.rev.1.bt2
reference_genome.rev.2.bt
```

This bowtie2 indexing procedure should be repeated with any genomes of contaminants to screen against.

## 2) Quality control (optional)

This section presents some optional steps to perform quality checks on the  `.fastq` reads sequence files. 

### <a id="fastqc">2.1) Quality report</a>

The [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) tool is used to produce quality reports on the read files. It will output a number of files in a target directory.

**Example fastqc command:**
```shell 
fastqc \
  -t $THREADS \ 
  -d tempDir/ \
  -o outputDir/ \ 
  reads.fastq
```
To be run for each `.fastq` reads file in the analysis.

### <a id="fastqscreen">2.2) Contamination screening</a>

The fastq_screen tool is used (with [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)) to detect and report the presence of contaminant sequences. A custom configuration file such as [this example](fastq_screen.conf) is necessary to use more than the default sequence files for detection.

Each of the `DATABASE` entries in the configuration file must be the path to a folder containing the index files produced by the same aligner tool as the `--aligner` argument indicates. In our case, we will use bowtie2. This is done as explained in the [reference genome indexing section](#bt2indexing). One of these `DATABASE` entries should be the index of the reference genome.

**Example fastq_screen command:**
```shell 
fastq_screen \
    --force \
    --subset 100000 \
    --threads $THREADS \
    --outdir outputDir/ \ 
    --aligner bowtie2 \
    --conf fastq_screen.conf \
    sample.fastq
```

To be run for each `.fastq` reads file in the analysis. 

## 3) Reads data processing and alignment

The pipeline for producing `.bam` alignment files from fastq sequencing files. For each map to produce, one fastq file is needed for single-end mapping while two files (one for the first mate and one for the second mate sequencing) are required for pair-end mapping. The input files may used in a compressed `.gz` archive.

### 3.1) Read trimming

The [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) tool is used to trim the reads. The command is different between single-end and pair-end applications.

**Example single-end trimmomatic command:**
```shell 
trimmomatic SE \
    -threads $THREADS \
    -phred33 \
    sample.fastq \ 
    sample.trimmed.fastq \
    ILLUMINACLIP:adapters.fasta:2:30:10 \
    LEADING:5 TRAILING:5 MINLEN:20
```

This tool requires the adapters sequences to be provided. The numbers trailing after the name of the adapters file are the defaults (see the [documentation](http://www.usadellab.org/cms/?page=trimmomatic) for details). The `LEADING` and `TRAILING` values are phred score threshold for trimming bases. MINLEN is the the length of reads after trimming under which the reads are discarded from the dataset.

**Example pair-end trimmomatic command:**
```shell 
trimmomatic PE \
    -threads $THREADS \
    -phred33 \
    -validatePairs \
    sample_1.fastq \
    sample_2.fastq \
    sample_1.trimmed.fastq \
    sample_1.trimmed.unpaired.fastq \
    sample_2.trimmed.fastq \
    sample_2.trimmed.unpaired.fastq \
    ILLUMINACLIP:adapters.fasta:2:30:10 \
    LEADING:5 TRAILING:5 MINLEN:20
```

The files marked unpaired are unnecessary and may be deleted immediately:

```shell 
rm -f sample_*.trimmed.unpaired.fastq
```

Optionally, the fastq files produced by trimmomatic may be compressed to save space: 

```shell 
pigz -p $THREADS -v sample.trimmed.fastq
```

To be run for each `.fastq` reads file in the analysis.

To control the quality of the data after trimming, the [fastqc command](#fastqc) may be repeated on the trimmed files.

### 3.2) Mapping

[STAR](https://github.com/alexdobin/STAR) is used to map reads to the genome. It allows aligning reads to the genome with the possibility of aligning the read over several sections separated by gaps which would be introns of the gene corresponding to the transcript from which that fragment was sequenced. The `--alignIntronMin` and `--alignIntronMax` arguments determine the range of lengths allowed for these gaps. 

The `--readFilesIn` argument is used to indicate the input read file. If two files are given the aligner will operate in pair-end mode. The `--readFilesCommand` argument should be adapted for the type of input (*e.g.* `cat` for an uncompressed `.fasta` file or `zcat` for a compressed `.fasta.bz` file). The `--genomeDir` argument should be the folder in which the [index](#indexing) was produced. The `--outSAMtype` argument determines the output type, in this case a sorted `.bam` file. The `--outFileNamePrefix` argument is a string that prefixes the output file names. The alignment file given the output options will be `<prefix>_Aligned.sortedByCoord.out.bam`.

The `--outFilterMismatchNmax` argument limits the number of allowed mismatches in each alignment. The `--outSAMmultNmax` determines how many alignments may be given for each read and `--outMultimapperOrder` determines how alignments are picked out of others with equal quality. Finally, the `--outFilterMultimapNmax` may be added to determine for how many alignments a multimapping read will be filtered out entirely. Note that this filter is applied by default with a value of 10 even if the argument is not specified. More information can be found in the [STAR manual](https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf).

Finally we index the alignment map with [samtools index](http://www.htslib.org/doc/samtools-index.html).

**Example single-end STAR alignment command:**
<!---
```shell
samtools index -@ $THREADS sample_Aligned.sortedByCoord.out.bam
```
--->

```shell
STAR --alignIntronMin 4 --alignIntronMax 16000 --runThreadN $THREADS \
  --readFilesIn sample.trimmed.fastq --readFilesCommand zcat \
  --genomeDir STAR_index/ \
  --outTmpDir $TMP \
  --outSAMtype BAM SortedByCoordinate --outFileNamePrefix sample \
  --outFilterMismatchNmax 2 --outSAMmultNmax 1 --outMultimapperOrder Random

sambamba index -t $THREADS sample_Aligned.sortedByCoord.out.bam
```

**Example paired-end STAR alignment command:**
<!---
```shell
samtools index -@ $THREADS sample_Aligned.sortedByCoord.out.bam
```
--->

```shell
STAR --alignIntronMin 4 --alignIntronMax 16000 --runThreadN $THREADS \
  --readFilesIn sample_1.trimmed.fastq sample_2.trimmed.fastq --readFilesCommand zcat \
  --genomeDir STAR_index/ \
  --outTmpDir $TMP \
  --outSAMtype BAM SortedByCoordinate --outFileNamePrefix sample \
  --outFilterMismatchNmax 2 --outSAMmultNmax 1 --outMultimapperOrder Random

sambamba index -t $THREADS sample_Aligned.sortedByCoord.out.bam
```


### 3.3) <a id="filtering">Reads filtering

<!---[samtools view](https://www.htslib.org/doc/samtools-view.html) --->
Filters are applied using [sambamba view](https://lomereiter.github.io/sambamba/docs/sambamba-view.html). The flag used to filter ecompasses the following:

|2572||
|:----|:--------------------------------------------|
| = 4 | read unmapped |
| + 8 | pair unmapped (does nothing in single-end) |
| + 512 | read fails platform/vendor quality checks |
| + 2048 | supplementary alignment  |

**Example filtering command:**
<!---
```shell
samtools view \
    -hb \
    -@ $THREADS \
    -F 2572 \
    -o sample.filtered.bam \
    sample_Aligned.sortedByCoord.out.bam
  
samtools index -@ $THREADS sample.filtered.bam
```
--->

```shell
sambamba view \
    -h -f bam \
    -t $THREADS \
    --num-filter /2572 \
    -o sample.filtered.bam \
    sample_Aligned.sortedByCoord.out.bam

sabamba index -t $THREADS sample.filtered.bam
```

To be run for each sample in the analysis. 

### 3.4) Masking ? 
TBD

### 3.6) Cleaning up (optional)

The pipeline, as presented, creates heavy intermediate files which allows for easy backtracking but is not economical in terms of disk space.

The only map files which are necessary for any further work are the files produced at [the filtering step](#filtering). All previous `.bam` files, their associated `.bai` index files and the trimmed `.fastq` files may be deleted or archived.

It is however advised to keep intermediate files if possible to be able to resume from any step if needed.

