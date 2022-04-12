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
    return data
            

# read cuts
def read_cuts(fr):
    data = {}
    for l in fr:
        ws = l.rstrip('\n').split('\t')
        chrom = ws[0]
        loc = int(ws[1])
        if chrom not in data:
            data[chrom] = []
        data[chrom].append(loc)
    return data


# calc 3' and 5' distance for a single entry
def add35(dist3, dist5, p, c, n, st):
    if st == '-':
        dist3.append(c-p)
        dist5.append(n-c)
    else:
        dist3.append(n-c)
        dist5.append(c-p)

# get distance
def calc_dist(rNMPs, res):
    dist3 = {'nuc':[], 'chrM':[]}
    dist5 = {'nuc':[], 'chrM':[]}
    for k, v in rNMPs.items():
        for st in ['+', '-']:
            name = 'chrM' if k == 'chrM' else 'nuc'
            data = sorted(v[st])
            i = 0
            j = 0
            cache = 0
            while i < len(data) and j < len(res[k]):
                if res[k][j] >= data[i]:
                    add35(dist3[name], dist5[name], cache, data[i], res[k][j], st)
                    i += 1
                else:
                    cache = res[k][j]
                    j += 1
    return dist3, dist5


# draw distance distribution
def draw(dists, out):
    sns.set(style='ticks')
    for k, v in dists.items():
        if not v:
            continue
        fig, ax = plt.subplots(figsize=(10,6))
        if k == 'nuc':
            sns.histplot(v, stat='probability', bins=50, binrange=(0, 200), ax=ax)
        else:
            sns.histplot(v, stat='probability', bins=50, binrange=(0, 200), ax=ax)
        sns.despine()
        ax.set_xlabel('Distance(nt)')
        plt.title(f'Distribution of distance to neareat RE cuts, {k}')
        fig.savefig(f'{out}_{k}.png')


def main():
    parser = argparse.ArgumentParser(description='Check distance between rNMPs and nearest RE site')
    parser.add_argument('bed', type=argparse.FileType('r'), help='rNMP incorporation bed file')
    parser.add_argument('re', type=argparse.FileType('r'), help='RE cut sites')
    parser.add_argument('-o', default='distance_distribution', help='Output base name')
    args = parser.parse_args()

    # get RE cuts and rNMPs
    rNMPs = read_bed(args.bed)
    res = read_cuts(args.re)

    # distances
    dist3, dist5 = calc_dist(rNMPs, res)

    # draw distance
    draw(dist3, args.o + '_3prime')
    draw(dist5, args.o + '_5prime')

    print('Done!')

if __name__ == '__main__':
    main()
