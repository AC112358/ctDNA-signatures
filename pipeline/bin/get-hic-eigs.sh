#!/bin/bash

set -e
if [ "$#" -ne 5 ]; then
    echo "usage: get-eigs.sh <assembly:hg38 or hg19> <fasta> <output> <regions> <coolfile>";
    exit 2;
fi
# preset args
assembly=$1
fasta=$2

# provided args
output=$3
regions=$4
cool_file=$5

output_res_kb=20
output_binsize=$((1000*$output_res_kb))
cool_handle="$cool_file::/resolutions/$output_binsize"

tempdir=$(mktemp -d)

echo "Writing view windows ..."
pythoncommand="import bioframe; import cooler; clr = cooler.Cooler('$cool_handle'); chromsizes=clr.chromsizes; cens=bioframe.fetch_centromeres('$assembly'); view=bioframe.make_chromarms(chromsizes, cens); view=view[view.chrom.astype('str').isin(clr.chromnames)].reset_index(drop=True).to_csv('$tempdir/view.tsv', index=False, header=False, sep='\t')"
python -c "$pythoncommand"

echo "Getting genome bins ..."
cooler dump --header -t bins $cool_handle | cut -f1-3 > $tempdir/bins.tsv

echo "Calculating GC content ..."
cooltools genome gc $tempdir/bins.tsv $fasta | awk -v OFS="\t" '{ if ( NF==4 ){ print $0 } }' > $tempdir/gc.tsv
bedtools intersect -a $tempdir/view.tsv -b $tempdir/gc.tsv -u -wa > $tempdir/filtered_view.tsv

echo "Calculating phased compartment eigs ..."
cooltools eigs-cis \
    -o $tempdir/output \
    --view $tempdir/filtered_view.tsv \
    --phasing-track $tempdir/gc.tsv \
    --n-eigs 1 \
    $cool_handle

awk -v OFS="\t" 'NF==5 && NR>1 {print $1,$2,$3,$5}' $tempdir/output.cis.vecs.tsv | sort -k1,1 -k2,2n > $output
rm -r $tempdir
