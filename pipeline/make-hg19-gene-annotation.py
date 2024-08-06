#!/bin/bash
set -e

gtf="gencode.v46lift37.basic.annotation.gff3.gz"
wget -O $gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh37_mapping/gencode.v46lift37.basic.annotation.gff3.gz
gzip -d $gtf.gz

query-gtf \
    --is-gff \
    -i $gtf \
    -type transcript \
    -attr tag \
    -vals MANE_Select \
    --contains \
    --header \
    -f "{chrom}\t{start}\t{end}\t{attributes[gene_name]}\t{attributes[gene_id]}\t{strand}\n"

rm $gtf