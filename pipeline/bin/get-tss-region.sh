#!/bin/bash

if [ $# -ne 4 ]; then
    echo -e "Invalid number of arguments. Expected 4 arguments.\nusage: get-tss-region.sh <promoter_len> <out> <regions> <genes_file>"
    exit 2
fi

promoter_len=$1
out=$2
regions=$3
genes=$4

set -e

awk -v OFS="\t" -v promlen=$promoter_len \
    '{ \
        ($6=="-") ? start = $3-promlen : start = $2-promlen; \
        ($6=="-") ? strand = "+" : strand = "-"; \
        (start < 0) ? start = 0 : start = start; \
    };
    {print $1,start,start+(2*promlen),"Promoter"}' $genes \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - -c 4 -o distinct \
    | sort -k1,1 -k2,2n \
    > $out
