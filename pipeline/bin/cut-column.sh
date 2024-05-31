#!/bin/bash

set -e
if [ $# -ne 4 ]; then
    echo -e "usage: cut-column.sh <column> <out> <regions> <bedfile>"
    exit 4
fi
col=$1
out=$2
bed=$4

cut -f "1-3,$col" $bed | sort -k1,1 -k2,2n > $out
