#!/bin/bash

if [ $# -ne 3 ]; then
    echo "Invalid number of arguments. Expected 3 arguments."
    exit 1
fi
out=$1
regions=$2
bed=$3

bedtools map -a $regions -b $bed \
    -c 5 -o mean -null nan -split | \
    awk '{print $NF}' > $out