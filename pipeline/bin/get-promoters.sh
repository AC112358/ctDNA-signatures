#!/bin/bash
if [ $# -ne 4 ]; then
    echo -e "Invalid number of arguments. Expected 4 arguments.\nusage: get-promoters.sh <promoter_len> <out> <regions> <genes>"
    exit 1
fi
promoter_len=$1
out=$2
regions=$3
genes=$4

set -e
promoter_regions=$(mktemp)

awk -v OFS="\t" -v promlen=$promoter_len \
    '{ \
        ($6=="-") ? start = $3 : start = $2-promlen; \
        ($6=="-") ? strand = "+" : strand = "-"; \
    };
    {print $1,start,start+promlen,$4,$5,strand}' $genes \
    | sort -k1,1 -k2,2n \
    | awk -v OFS="\t" '$2>=0' \
    | bedtools subtract -a - -b $genes -sorted \
    | sort -k1,1 -k2,2n \
    > $promoter_regions

bedtools subtract \
    -a $promoter_regions \
    -b $promoter_regions \
    -A -f 0.25 -r -e -S -sorted \
    | bedtools merge -i - -s -c 6 -o distinct > $out
