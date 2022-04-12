#!/usr/bin/env python3

import argparse
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict


# read bed file
def read_bed(fr):
    data = {}
    for l in fr:
        ws = l.rstrip('\n').split('\t')
        chrom = ws[0]
        st = ws[5]
        loc = int(ws[2])
        if chrom not in data:
            data[chrom] = {'+':[], '-':[]}
        data[chrom][st].append(loc)
    dists = {'nuc':[], 'chrM':[]}
    for k, v in data.items():
        for st in ['+', '-']:
            v[st].sort()
            for i in range(len(v[st])-1):
                d = v[st][i+1] - v[st][i]
                # skip rNMPs at the same location
                if not d:
                    continue
                # nuc and chrM separate
                if k == 'chrM':
                    dists[k].append(d)
                else:
                    dists['nuc'].append(d)
    return dists
            

# draw distance distribution
def draw(dists, out):
    sns.set(style='ticks')
    for k, v in dists.items():
        if not v:
            continue
        fig, ax = plt.subplots(figsize=(10,6))
        if k == 'nuc':
            sns.histplot(v, stat='probability', bins=50, binrange=(0, 50), ax=ax)
        else:
            sns.histplot(v, stat='probability', bins=50, binrange=(0, 50), ax=ax)
        sns.despine()
        ax.set_xlabel('Distance(nt)')
        plt.title(f'Distribution of gap length between two rNMPs, {k}')
        fig.savefig(f'{out}_{k}.png')


def main():
    parser = argparse.ArgumentParser(description='Check distance between rNMPs')
    parser.add_argument('bed', type=argparse.FileType('r'), help='rNMP incorporation bed file')
    parser.add_argument('-o', default='distance_distribution', help='Output base name')
    args = parser.parse_args()

    # get distances between rNMPs
    dists = read_bed(args.bed)

    # draw distance
    draw(dists, args.o)

    print('Done!')

if __name__ == '__main__':
    main()
