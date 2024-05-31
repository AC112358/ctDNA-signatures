#!/bin/bash
set -e

gtf="MANE.GRCh38.v1.3.ensembl_genomic.gtf"
wget -O $gtf.gz https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.3/MANE.GRCh38.v1.3.ensembl_genomic.gtf.gz
gzip -d $gtf.gz

python bin/query_gtf.py \
    -i $gtf \
    -type transcript \
    -attr gene_id \
    --header \
    -f "{chrom}\t{start}\t{end}\t{attributes[gene_name]}\t{attributes[gene_id]}\t{strand}\n"

rm $gtf