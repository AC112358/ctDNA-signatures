#!/usr/bin/env python3
import pandas as pd
import sys
import argparse

def main(output=sys.stdout,*,
        encode_quantification_file,
        curated_genes_file,
        ):
    
    # read and reformate the quantification file
    gene_annotation = pd.read_csv(curated_genes_file, sep = '\t')
    gene_annotation.columns = [x.strip('#').removeprefix('attributes.') for x in gene_annotation.columns]
    gene_annotation.gene_id = gene_annotation.gene_id.str.split('.').str[0]

    quantifications = pd.read_csv(encode_quantification_file, sep = '\t')
    quantifications.gene_id = quantifications.gene_id.str.split('.').str[0]
    #quantifications.transcript_id = quantifications.transcript_id.str.split('.').str[0]

    gene_annotation = gene_annotation.merge(
        quantifications[['gene_id','pme_TPM']], 
        on = 'gene_id',
        how = 'inner'
    )

    print(f'Found {len(gene_annotation)} shared genes...', file=sys.stderr)

    #expression_by_id = gene_annotation.groupby('gene_id')['pme_FPKM'].sum()
    # select the most expressed transcript for each gene
    #gene_annotation = gene_annotation.sort_values('pme_FPKM', ascending=False)\
    #                    .groupby('gene_id').first() 
    
    #gene_annotation = gene_annotation.join(expression_by_id, on='gene_id', rsuffix='_sum').reset_index()
    
    bed_features = ['chrom','start','end']
    gene_annotation.sort_values(bed_features, inplace=True)
    
    # save the gene expression
    gene_annotation[[*bed_features, 'gene_id','pme_TPM', 'strand']]\
        .to_csv(output, sep = '\t', index = None, header = None)


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--encode-quantification-file','-quant', type=str, required=True)
    parser.add_argument('--curated-genes-file','-genes', type=str, required=True)
    parser.add_argument('--output','-o', type=argparse.FileType('w'), default=sys.stdout)
    args = parser.parse_args()
    main(**vars(args))
