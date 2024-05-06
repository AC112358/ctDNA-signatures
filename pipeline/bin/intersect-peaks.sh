#!/bin/bash

if [ $# -lt 4 ]; then
    echo -e "Invalid number of arguments. Expected at least 4 arguments."
    exit 1
fi
feature_name=$1
out=$2
regions=$3
bedfiles=${@:4}

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
