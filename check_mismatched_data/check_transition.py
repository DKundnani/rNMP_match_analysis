#!/usr/bin/env python3

import argparse
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import pandas as pd


# read bed file
def read_fa(fr, f):
    res = {x:[0]*(2*f) for x in ['Same', 'Different']}
    for l in fr:
        if l[0] == '>':
            continue
        l = l.upper()
        for i in range(2*f):
            if l[i] == 'N' or l[i+1] == 'N':
                continue
            if l[i] == l[i+1]:
                res['Same'][i] += 1
            else:
                res['Different'][i] += 1
    res = pd.DataFrame(res)
    res['Total'] = res.sum(axis=1)
    res['Location'] = list(range(-f, f))
    return res
            

# draw distance distribution
def draw(df, out, f):
    sns.set(style='ticks', font_scale=2)
    fig, ax = plt.subplots(figsize=(12,5))
    plt.subplots_adjust(bottom=.15, left=.1, right=0.8)
    pal = ['r', 'b']
    s = [0] * len(df)
    for i, c in enumerate(['Same', 'Different']):
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
    parser.add_argument('-o', default='transition', help='Output sequences')
    args = parser.parse_args()

    # get composition at each location
    freqs = read_fa(args.fa, args.f)

    # draw distance
    draw(freqs, args.o, args.f)

    print('Done!')

if __name__ == '__main__':
    main()
