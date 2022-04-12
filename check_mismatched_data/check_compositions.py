#!/usr/bin/env python3

import argparse
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import pandas as pd
import itertools as it
import numpy as np

# generate possible combinations
def generate_combination(k):
    for x in it.product(list('ACGT'), repeat=k):
        yield ''.join(x)


# read bed file
def read_fa(fr, f, k):
    res = {x:[0]*(2*f+2-k) for x in generate_combination(k)}
    for l in fr:
        if l[0] == '>':
            continue
        for i in range(2*f+2-k):
            c = l[i:i+k]
            if c in res:
                res[c][i] += 1
    res = pd.DataFrame(res)
    res['Total'] = res.sum(axis=1)
    res['Location'] = list(range(-f, f+2-k))
    return res
            

# draw distance distribution
def draw(df, out, f, k):
    sns.set(style='ticks', font_scale=2)
    fig, ax = plt.subplots(figsize=(12,5))
    plt.subplots_adjust(bottom=.15, left=.1, right=0.8)
    if k == 1:
        pal = ['#E31A1C', '#1F78B4', '#FFFFB9', '#33A02C']
    else:
        cache = list(sns.hls_palette(len(df.columns)-2))
        pal = []
        for i in range(len(cache)//2):
            pal.append(cache[i])
            pal.append(cache[len(cache)//2+i])
    s = [0] * (len(df))
    for i, c in enumerate(generate_combination(k)):
        ax.bar(df.Location, df[c]/df['Total'], label=c, \
            color=pal[i], bottom=s, linewidth=0, width=1)
        s = df[c]/df['Total']+s
    plt.legend(bbox_to_anchor=(1,1))
    sns.despine()
    ax.set_xlabel('Location')
    ax.set_ylabel('Fraction')
    ax.set_xlim((-f,f))
    fig.savefig(f'{out}.png')


def main():
    parser = argparse.ArgumentParser(description='Check composition at each location')
    parser.add_argument('fa', type=argparse.FileType('r'), help='Fasta sequences')
    parser.add_argument('-f', type=int, default=20, help='Flank length (20 nt)')
    parser.add_argument('-k', type=int, default=1, help='K-mer length, default=1')
    parser.add_argument('-o', default='composition', help='Output sequences')
    args = parser.parse_args()

    # get composition at each location
    freqs = read_fa(args.fa, args.f, args.k)

    # draw distance
    draw(freqs, args.o, args.f, args.k)

    print('Done!')

if __name__ == '__main__':
    main()
