#!/bin/bash

if [ $# -lt 5 ]; then
    echo -e "Invalid number of arguments. Expected at least 5 arguments."
    echo -e "usage: intersect-peaks.sh <feature_name> <threshold> <out> <regions> <bedfiles>"
    exit 2
fi
feature_name=$1
threshold=$2
out=$3
regions=$4
bedfiles=${@:5}

set -e
tmpdir=$(mktemp -d)

for file in $bedfiles; do 
    bedtools sort -i $file > "$tmpdir/$(basename $file).sorted.bed"; 
done

bedtools multiinter -i $tmpdir/*.sorted.bed \
    | bedtools sort -i - \
    | bedtools merge -i - -c 4 -o max \
    | awk -v OFS="\t" -v threshold=2 -v featurename=$feature_name '$4>=threshold {print $1,$2,$3,featurename}' \
    > $out
