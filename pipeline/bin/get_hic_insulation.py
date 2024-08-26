#!/usr/bin/env python

import sys
import bioframe
import cooler
from cooltools import insulation

def main(cmd, coolfile):
    resolution=20000
    cool_handle=f'{coolfile}::/resolutions/{resolution}'
    cooler_obj=cooler.Cooler(cool_handle)
    windows = [resolution]

    insulation_table = insulation(cooler_obj, windows, verbose=True)

    insulation_table = insulation_table[(~insulation_table['is_bad_bin']) & (insulation_table.chrom != 'MT')]
    insulation_table['chrom'] = insulation_table.chrom.astype(str)
    
    if cmd=='score':
        col='log2_insulation_score_' + str(resolution)
        print(
            insulation_table.head(10)
        )
        insulation_table = insulation_table[~insulation_table[col].isna()][['chrom','start','end',col]]

    elif cmd=='signpost':
        insulation_table = insulation_table[insulation_table['is_boundary_' + str(resolution)]][['chrom','start','end']]
        insulation_table['boundary']='Boundary'
    else:
        raise ValueError(f'Command {cmd} is not an option!')
    
    insulation_table = insulation_table.sort_values(['chrom','start'])
    
    return insulation_table
        

if __name__=="__main__":

    args=sys.argv[1:]
    if not len(args)==4:
        print('Usage: insulation.py <command: score or signpost> <output> <regions> <input>', file=sys.stderr)
        sys.exit(2)
    
    main(args[0], args[-1]).to_csv(args[1],sep='\t', index=False, header=False)