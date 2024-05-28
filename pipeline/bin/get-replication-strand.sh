#!/bin/bash

set -e
if [ $# -lt 2 ]; then
    echo -e "Invalid number of arguments. Expected at least 4 arguments."
    exit 3
fi
out=$1
chromsizes=$2
bigwigs=${@:2}

smooth_bigwig=$(mktemp)
windows=$(mktemp)
binned_counts=$(mktemp)

./average-bigwigAverage $smooth_bigwig $bigwigs

sort -k1,1 $chromsizes \
    | bedtools makewindows -g - -w 100000 -s 25000 -i srcwinnum > $windows

bigWigAverageOverBed -bedOut=$binned_counts.bed $smooth_bigwig $windows $binned_counts
rm $binned_counts

python replication_hmm.py $binned_counts.bed > $out
