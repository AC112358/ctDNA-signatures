#!/usr/bin/env python3
import pandas as pd
import sys
import argparse

def main(output=sys.stdout,*,
        encode_quantification_file,
        curated_genes_file,
        ):
    
    # read and reformate the quantification file
    quantifications = pd.read_csv(encode_quantification_file, sep = '\t')

    gene_annotation = pd.read_csv(curated_genes_file, sep = '\t')
    gene_annotation.columns = gene_annotation.columns.str.split('.').str[-1]

    gene_annotation = gene_annotation.merge(
        quantifications[['transcript_id','gene_id','pme_FPKM']], 
        left_on = 'value',
        right_on = 'transcript_id',
        how = 'inner'
    )

    # select the most expressed transcript for each gene
    gene_annotation = gene_annotation.sort_values('pme_FPKM', ascending=False)\
                        .groupby('gene_id').first().sort_values(['chrom','txStart']).reset_index()
    
    gene_annotation.sort_values(['chrom','txStart'], inplace=True)

    bed_features = ['chrom','txStart','txEnd']

    # save the gene expression
    gene_annotation[[*bed_features, 'geneSymbol','pme_FPKM', 'strand']]\
        .to_csv(output, sep = '\t', index = None, header = None)


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--encode-quantification-file','-quant', type=str, required=True)
    parser.add_argument('--curated-genes-file','-genes', type=str, required=True)
    parser.add_argument('--output','-o', type=argparse.FileType('w'), default=sys.stdout)
    args = parser.parse_args()
    main(**vars(args))

    