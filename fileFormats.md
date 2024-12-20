# File formats quick reference

A reference page of the different file formats used by genomics tools mentionned in the pipeline guides of this repository.

## FASTA
**Extensions:** `.fasta`, `.fa`

**Type:** Text.

**Description:** A text-based format for representing either nucleotide sequences or amino acid (protein) sequences, in which nucleotides or amino acids are represented using single-letter codes. A sequence begins with a greater-than character (">") followed by a description of the sequence (all in a single line). The lines immediately following the description line are the sequence representation. A sequence may be split between any number of lines.

**Example:**
```
>AT1G01010.1
GTCTTCCTCCCTCCAAATTATTAGATATACCAAACCAGAGAAAACAAATACATAATCGGAGAAATACAGATTACAGAGAGCGAGAGAGATCGACGGCGAAGCTCTTTACCCGGAAACCATTGAAATCGGACGGTTTAGTGAAAATGGAGGATCAAGTTGGGTTTGGGTTCCGTCCGAACGACGAGGAGCTCGTTGGTCACTATCTCCGTAACAAAATCGAAGGAAACACTAGCCGCGACGTTGAAGTAGCCATCAGCGAGGTCAACATCTGTAGCTACGATCCTTGGAACTTGCGCTTCCAGTCAAAGTACAAATCGAGAGATGCTATGTGGTACTTCTTCTCTCGTAGAGAAAACAACAAAGGGAATCGACAGAGCAGGACAACGGTTTCTGGTAAATGGAAGCTTACCGGAGAATCTGTTGAGGTCAAGGACCAGTGGGGATTTTGTAGTGAGGGCTTTCGTGGTAAGATTGGTCATAAAAGGGTTTTGGTGTTCCTCGATGGAAGATACCCTGACAAAACCAAATCTGATTGGGTTATCCACGAGTTCCACTACGACCTCTTACCAGAACATCAGAGGACATATGTCATCTGCAGACTTGAGTACAAGGGTGATGATGCGGACATTCTATCTGCTTATGCAATAGATCCCACTCCCGCTTTTGTCCCCAATATGACTAGTAGTGCAGGTTCTGTGGTGAGTCTTTCTCCATATACACTTAGCTTTGAGTAGGCAGATCAAAAAAGAGCTTGTGTCTACTGATTTGATGTTTTCCTAAACTGTTGATTCGTTTCAGGTCAACCAATCACGTCAACGAAATTCAGGATCTTACAACACTTACTCTGAGTATGATTCAGCAAATCATGGCCAGCAGTTTAATGAAAACTCTAACATTATGCAGCAGCAACCACTTCAAGGATCATTCAACCCTCTCCTTGAGTATGATTTTGCAAATCACGGCGGTCAGTGGCTGAGTGACTATATCGACCTGCAACAGCAAGTTCCTTACTTGGCACCTTATGAAAATGAGTCGGAGATGATTTGGAAGCATGTGATTGAAGAAAATTTTGAGTTTTTGGTAGATGAAAGGACATCTATGCAACAGCATTACAGTGATCACCGGCCCAAAAAACCTGTGTCTGGGGTTTTGCCTGATGATAGCAGTGATACTGAAACTGGATCAATGATTTTCGAAGACACTTCGAGCTCCACTGATAGTGTTGGTAGTTCAGATGAACCGGGCCATACTCGTATAGATGATATTCCATCATTGAACATTATTGAGCCTTTGCACAATTATAAGGCACAAGAGCAACCAAAGCAGCAGAGCAAAGAAAAGGTGATAAGTTCGCAGAAAAGCGAATGCGAGTGGAAAATGGCTGAAGACTCGATCAAGATACCTCCATCCACCAACACGGTGAAGCAGAGCTGGATTGTTTTGGAGAATGCACAGTGGAACTATCTCAAGAACATGATCATTGGTGTCTTGTTGTTCATCTCCGTCATTAGTTGGATCATTCTTGTTGGTTAAGAGGTCAAATCGGATTCTTGCTCAAAATTTGTATTTCTTAGAATGTGTGTTTTTTTTTGTTTTTTTTTCTTTGCTCTGTTTTCTCGCTCCGGAAAAGTTTGAAGTTATATTTTATTAGTATGTAAAGAAGAGAAAAAGGGGGAAAGAAGAGAGAAGAAAAATGCAGAAAATCATATATATGAATTGGAAAAAAGTATATGTAATAATAATTAGTGC
>AT1G01010.2
AAATTATTAGATATACCAAACCAGAGAAAACAAATACATAATCGGAGAAATACAGATTACAGAGAGCGAGAGAGATCGACGGCGAAGCTCTTTACCCGGAAACCATTGAAATCGGACGGTTTAGTGAAAATGGAGGATCAAGTTGGGTTTGGGTTCCGTCCGAACGACGAGGAGCTCGTTGGTCACTATCTCCGTAACAAAATCGAAGGAAACACTAGCCGCGACGTTGAAGTAGCCATCAGCGAGGTCAACATCTGTAGCTACGATCCTTGGAACTTGCGCTTCCAGTCAAAGTACAAATCGAGAGATGCTATGTGGTACTTCTTCTCTCGTAGAGAAAACAACAAAGGGAATCGACAGAGCAGGACAACGGTTTCTGGTAAATGGAAGCTTACCGGAGAATCTGTTGAGGTCAAGGACCAGTGGGGATTTTGTAGTGAGGGCTTTCGTGGTAAGATTGGTCATAAAAGGGTTTTGGTGTTCCTCGATGGAAGATACCCTGACAAAACCAAATCTGATTGGGTTATCCACGAGTTCCACTACGACCTCTTACCAGAACATCAGAGGACATATGTCATCTGCAGACTTGAGTACAAGGGTGATGATGCGGACATTCTATCTGCTTATGCAATAGATCCCACTCCCGCTTTTGTCCCCAATATGACTAGTAGTGCAGGTTCTGTGGTCAACCAATCACGTCAACGAAATTCAGGATCTTACAACACTTACTCTGAGTATGATTCAGCAAATCATGGCCAGCAGTTTAATGAAAACTCTAACATTATGCAGCAGCAACCACTTCAAGGATCATTCAACCCTCTCCTTGAGTATGATTTTGCAAATCACGGCGGTCAGTGGCTGAGTGACTATATCGACCTGCAACAGCAAGTTCCTTACTTGGCACCTTATGAAAATGAGTCGGAGATGATTTGGAAGCATGTGATTGAAGAAAATTTTGAGTTTTTGGTAGATGAAAGGACATCTATGCAACAGCATTACAGTGATCACCGGCCCAAAAAACCTGTGTCTGGGGTTTTGCCTGATGATAGCAGTGATACTGAAACTGGATCAATGATTTTCGAAGACACTTCGAGCTCCACTGATAGTGTTGGTAGTTCAGATGAACCGGGCCATACTCGTATAGATGATATTCCATCATTGAACATTATTGAGCCTTTGCACAATTATAAGGCACAAGAGCAACCAAAGCAGCAGAGCAAAGAAAAGGTGATAAGTTCGCAGAAAAGCGAATGCGAGTGGAAAATGGCTGAAGACTCGATCAAGATACCTCCATCCACCAACACGGTGAAGCAGAGCTGGATTGTTTTGGAGAATGCACAGTGGAACTATCTCAAGAACATGATCATTGGTGTCTTGTTGTTCATCTCCGTCATTAGTTGGATCATTCTTGTTGGTTAAGAGGTCAAATCGGATTCTTGCTCAAAATTTGTATTTCTTAGAATGTGTGTTTTTTTTTGTTTTTTTTTCTTTGCTCTGTTTTCTCGCTCCGGAAAAGTTTGAAGTTATATTTTATTAGTATGTAAAGAAGAGAAAAAGGGGGAAAGAAGAGAGAAGAAAAATGCAGAAAATCATATATATGAATTGGAAAAAAGTATATGTAATAATAATTAGTGCATCGTTTTGTGGTGTAGTTTATATAAATAAAGTGATATATAGTCTTGT
```

**Note:** Often compressed as a GZ archive file with the extension `fasta.gz` or `fa.gz`.

## FASTQ
**Extensions:** `.fastq`, `.fq`

**Type:** Text.

**Description:** A text-based format for storing both a biological sequence (usually nucleotide sequence) and its corresponding quality scores. Both the sequence letter and quality score are each encoded with a single ASCII character.

**Example:**
```
@V350134665L3C001R00100001579:0:0:0:0 1:N:0:CTCTAC
AGTCGAATATGACTTGATCTCATGTGTATGATTAAGTATAAGAACTTAAACCGCAACAAGATCTTAAAGGCGTAAGAATTGTATCCTTGTTAAAAGACAC
+
IIIIIIIIIIIIIIIIHIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFIII
@V350134665L3C001R00100001642:0:0:0:0 1:N:0:CTCTAC
AAATCACTCTTGACAGTAATATCTGTTGTATATGTAAATCCTAGATATGACAATATGCGGAATTCATCCATGAAAATGAAAAACAAGGGGGGGTTGCTAA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIEIIIII
```

**Note:** Often compressed as a GZ archive file with the extension `fastq.gz` or `fq.gz`.

## SAM/BAM
**Extensions:** `.sam`, `.bam`

**Type:** Text, TSV (tab-separated values) (SAM) or Binary (BAM). 

**Description:** Sequence Alignment Map (SAM) is a text-based format originally for storing biological sequences aligned to a reference sequence. Consists of a header and an alignment section. Alignment sections have 11 mandatory fields, as well as a variable number of optional fields.
```
Col 	Field 	Type 	Brief description
1       QNAME 	String 	Query template NAME
2       FLAG 	Int 	bitwise FLAG
3       RNAME 	String 	References sequence NAME
4       POS 	Int 	1- based leftmost mapping POSition
5       MAPQ 	Int 	MAPping Quality
6       CIGAR 	String 	CIGAR string
7       RNEXT 	String 	Ref. name of the mate/next read
8       PNEXT 	Int 	Position of the mate/next read
9       TLEN 	Int 	observed Template LENgth
10      SEQ 	String 	segment SEQuence
11      QUAL 	String 	ASCII of Phred-scaled base QUALity+33
```
A Binary Alignment Map (BAM) file stores the same data in a compressed binary representation.

**Full specification:** [SAMv1.pdf](https://samtools.github.io/hts-specs/SAMv1.pdf)

**Example (SAM):**
```
@HD VN:1.6 SO:coordinate
@SQ SN:ref LN:45
r001 99 ref 7 30 8M2I4M1D3M = 37 39 TTAGATAAAGGATACTG *
r002 0 ref 9 30 3S6M1P1I4M * 0 0 AAAAGATAAGGATA *
r003 0 ref 9 30 5S6M * 0 0 GCCTAAGCTAA * SA:Z:ref,29,-,6H5M,17,0;
r004 0 ref 16 30 6M14N5M * 0 0 ATAGCTTCAGC *
r003 2064 ref 29 17 6H5M * 0 0 TAGGC * SA:Z:ref,9,+,5S6M,30,1;
r001 147 ref 37 30 9M = 7 -39 CAGCGGCAT * NM:i:1
```

## SAI/BAI
**Extensions:** `.sam.sai`, `.bam.bai`

**Type:** Binary.

**Description:** Index file of a correspondig SAM or BAM file respectively.


## GFF/GTF
**Extensions:** `.gff`, `.gtf`

**Type:** Text, TSV (tab-separated values).

**Description:** The GFF (General Feature Format) format consists of one line per feature, each containing 9 columns of data, plus optional track definition lines. The GTF (General Transfer Format) is identical to GFF version 2.
  * **seqname** - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
  * **source** - name of the program that generated this feature, or the data source (database or project name)
  * **feature** - feature type name, e.g. Gene, Variation, Similarity
  * **start** - Start position* of the feature, with sequence numbering starting at 1.
  * **end** - End position* of the feature, with sequence numbering starting at 1.
  * **score** - A floating point value.
  * **strand** - defined as + (forward) or - (reverse).
  * **frame** - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
  * **attribute** - A semicolon-separated list of tag-value pairs, providing additional information about each feature.

**Full specification:** [gff3.md](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)

**Example:**
```
1 transcribed_unprocessed_pseudogene  gene        11869 14409 . + . gene_id "ENSG00000223972"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; 
1 processed_transcript                transcript  11869 14409 . + . gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-002"; transcript_source "havana";
```

## BED
**Extensions:** `.bed`

**Type:** Text, whitespace-delimited, preferably TSV (tab-separated values).

**Description:**The BED (Browser Extensible Data) format is a text file format used to store genomic regions as coordinates and associated annotations. A BED file must contain at least the first 3 fields of the specification table below A BED*n* file is a file which contains the first *n* specified fields. A BED*n*+*m* is a BED file starting with the first *n* specified fields followed by *m* custom fields. 

Track lines (the word "track" followed by a series of space separated key=value pairs) may be added at the start of the file to add display information relevant to a particular sequence browser but this is no longer technically a BED file.

<table align="left" width="1000" cellspacing="0" cellpadding="0">
    <tr>
        <td><img src="./images/bed_fields.png" align="left" width="1000px"></td>
    </tr>
</table>

**Full specification:** [BEDv1.pdf](https://samtools.github.io/hts-specs/BEDv1.pdf)

**Example:**
```
chr7 127471196 127472363 Pos1 0 + 127471196 127472363 255,0,0
chr7 127472363 127473530 Pos2 0 + 127472363 127473530 255,0,0
chr7 127473530 127474697 Pos3 0 + 127473530 127474697 255,0,0
chr7 127474697 127475864 Pos4 0 + 127474697 127475864 255,0,0
chr7 127475864 127477031 Neg1 0 - 127475864 127477031 0,0,255
```

## BedGraph
**Extensions:** `.bedgraph`

**Type:** Text, whitespace-delimited, preferably TSV (tab-separated values).

**Description:** A BedGraph file is similar to a BED3+1 file where the $4^{th}$ field is a score or other quantitative value. It also must start with a track lines containing `type=bedGraph`.

**Example:**
```
track type=bedGraph name="BedGraph Format" description="BedGraph format" priority=20
chr19 59302000 59302300 -1.0
chr19 59302300 59302600 -0.75
chr19 59302600 59302900 -0.50
chr19 59302900 59303200 -0.25
```

## BigWig
**Extensions:** `.bigwig`

**Type:** Binary.

**Description:** BigWig files are a compressed, indexed, binary format for genome-wide signal data for calculations (e.g. GC percent) or experiments (e.g. ChIP-seq/RNA-seq read depth).

