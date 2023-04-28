#!/bin/bash

##################################################################################################
#~ mergeOverlappingRegions                                                                      ~#
#~ Author: Adrien Vidal                                                                         ~#
#~ Last Edited: 10 mar. 2023                                                                    ~#
##################################################################################################

#~ Options:

function usage {
  echo "Summary: Searches for overlaps of the genomic regions of the between two files and generates merged regions bed file out of selected regions."
  echo "Usage: $(basename $0) -a <bed/gff/vcf/bam> -b <bed/gff/vcf/bam> [-l <min overlap>]"
  echo $'\t'"-a"$'\t'"file A. Any format accepted by bedtools. Each feature in A is compared to B in search of overlaps."
  echo $'\t'"-b"$'\t'"files B. Any format accepted by bedtools. Use -b mutiple times to search overlaps of a in more than one file."
  echo $'\t'"-l"$'\t'"overlap. Value between 0 and 1. Minimum overlap between two regions from A and B respectively as a fraction of either region."
}

while getopts ":a:b:l:" option
do
   case $option in
    a)
        fileA="$OPTARG";;
    b)
        fileB="$OPTARG";;
    l)
        overlap="$OPTARG";;
    \?) # Invalid option
        usage
        exit;;
   esac
done

#~ Tests:

if [[ -z "$fileA" ]]
then
    usage
    exit;
fi

if [[ -z "$fileB" ]]
then
    usage
    exit;
fi

if [[ -z $overlap ]]
then
    overlap=1E-9
fi

#~ Execution:

# Logic:
#
# file A -> chrA, startA, endA \                                                            | chrA, startA, endA |
#                               >--> | chrA, startA, endA, chrB, startB, endB, overlap | -> | ...                | -> sort -> merge
# file B -> chrB, startB, endB /     | ...                                             |    | chrB, startB, endB |

bedtools intersect -sorted -wa -wb -f $overlap -F $overlap -e \
    -a <(cut -f 1,2,3 $fileA | bedtools sort) \
    -b <(cut -f 1,2,3 $fileB | bedtools sort) \
| awk ' { print $1 "\t" $2 "\t" $3 }
        { print $4 "\t" $5 "\t" $6 }' \
| bedtools sort | bedtools merge