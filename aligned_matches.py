#!usr/bin/env python3
#Author: Deepali L. Kundnani
#Description:

import argparse
import os
import pandas as pd  

#fraction of matches out of total aligned
def fraction(m):
    if len(m)>0:
        return(m.count(1)/len(m))
    else:
        return('NA')
    
#Specifies Match or Mismatch    
def status(row):
    if row['seq_fa'][49] != row['seq_ref']:
        return("Mismatch")
    else:
        return("Match")

#Gets a list of aligned match(1) or mismatch(0) score within a distance of 50bp
def alignment(c,nc,f,col):
    dist = f.iloc[nc,:]['start'] - f.iloc[c,:]['start']
    match=[]
    while dist<50:
        if dist==0:
            nc=nc+1
            dist = f.iloc[nc,:]['start'] - f.iloc[c,:]['start']
            continue
    
        if f.iloc[c,:]['chr'] == f.iloc[nc,:]['chr'] and  f.iloc[c,:]['strand'] == f.iloc[nc,:]['strand']:
            if f.iloc[c,:]['seq_fa'][49]==f.iloc[nc,:]['seq_fa'][49-dist]:
                match.append(1)
            else:
                match.append(0)
        else:
            break
        
        nc = nc + 1
        if nc==(col):
            break
        else:
            dist = f.iloc[nc,:]['start'] - f.iloc[c,:]['start']    
    
    return(match)
    
    
def main():
    parser = argparse.ArgumentParser(description='example: python3 aligned_matches.py -s ~/scratch/stats/Infections -m ~/scratch/merge/mergedforPRS.vcf -o ~/scratch/stats/PRS -r ~/scratch/stats/percentile -c INFO_y')
    parser.add_argument('-f','--file', required=True ,help='Match Analysis output with upto 50 bp of raw fastq read sequence')
    parser.add_argument('-o','--out', default="aligned_matches_", help='Output folder/path(full path)')
    args = parser.parse_args()

    
    bd=pd.read_csv(args.file, sep='\t', header=None)
    bd.columns = ['chr', 'start', 'stop', 'seq_fa', 'seq_ref', 'strand']
    bd=bd.sort_values(['chr', 'strand', 'start'])

    o=open(args.out+args.file, 'w')
    for col in range(len(bd.index)-2):
        row=list(bd.iloc[col,:])
        row.append(status(bd.iloc[col,:]))
        match=alignment(col, col+1,bd,len(bd.index))
        row.append(len(match))
        row.append(fraction(match))
        o.write('\t'.join(str(x) for x in row))
        o.write('\n')
    o.close()    
     
if __name__ == '__main__':
    main() 