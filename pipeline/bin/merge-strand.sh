#!/bin/bash

if [ $# -ne 3 ]; then
    echo "Invalid number of arguments. Expected 3 arguments."
    exit 1
fi
out=$1
regions=$2
bed=$3

set -e
bedtools merge -i $bed -s -c 6 -o distinct > $out
