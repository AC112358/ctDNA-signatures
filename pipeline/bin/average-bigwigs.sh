#!/bin/bash

set -e
if [ $# -lt 2 ]; then
    echo -e "Invalid number of arguments. Expected at least 4 arguments."
    exit 3
fi
out=$1
bigwigs=${@:2}

bigwigAverage -b $bigwigs -o $out
