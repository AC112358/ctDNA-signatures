#!/usr/env/bin python3

import pandas as pd

genes=pd.read_csv('~/Downloads/hg19.knownCanonical.gtf', sep='\t')
genes.columns='attributes.gene_name,chrom,strand,start,end,attributes.gene_id'.split(',')
genes.drop_duplicates(subset='attributes.gene_id', inplace=True)
genes[['chrom','start','end','attributes.gene_name','attributes.gene_id','strand']]\
    .to_csv('hg19.gene.annotation.bed', sep='\t', index=None)
