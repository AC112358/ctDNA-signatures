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
    quantifications['ensemble_transcript_id'] = quantifications['transcript_id(s)'].str.split('.').str[0]

    gene_annotation = pd.read_csv(curated_genes_file, sep = '\t')

    gene_annotation = gene_annotation.merge(
        quantifications[['ensemble_transcript_id','pme_FPKM']], 
        left_on = 'hg19.knownToEnsembl.value',
        right_on = 'ensemble_transcript_id',
        how = 'inner'
    )
    gene_annotation.columns = gene_annotation.columns.str.split('.').str[-1]

    # select the most expressed transcript for each gene
    gene_annotation['gene_length'] = gene_annotation.txEnd - gene_annotation.txStart
    gene_annotation = gene_annotation.sort_values('pme_FPKM', ascending=False)\
                        .groupby('ensemble_transcript_id').first().sort_values(['chrom','txStart']).reset_index()
    
    gene_annotation.sort_values(['chrom','txStart'], inplace=True)

    bed_features = ['chrom','txStart','txEnd']

    # save the gene expression
    gene_annotation[[*bed_features, 'ensemble_transcript_id','pme_FPKM', 'strand']]\
        .to_csv(output, sep = '\t', index = None, header = None)
    

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--encode-quantification-file','-quant', type=str, required=True)
    parser.add_argument('--curated-genes-file','-genes', type=str, required=True)
    parser.add_argument('--output','-o', type=argparse.FileType('w'), default=sys.stdout)
    args = parser.parse_args()
    main(**vars(args))

    