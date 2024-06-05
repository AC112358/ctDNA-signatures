#!/bin/bash

set -e
if [ $# -lt 3 ]; then
    echo -e "Invalid number of arguments. Expected at least 3 arguments."
    exit 3
fi
out=$1
bigwigs=${@:3}

bigwigAverage -b $bigwigs -o $out
