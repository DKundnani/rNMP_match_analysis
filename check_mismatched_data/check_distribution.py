#!/usr/bin/env python3

import argparse
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict


# generate binsizes
def generate_bs(b, mb):
    data = defaultdict(lambda: b)
    data['chrM'] = mb
    return data


# read fasta index files
def read_fai(fr):
    data = {}
    for l in fr:
        ws = l.rstrip('\n').split('\t')
        data[ws[0]] = int(ws[1])
    return data


# load bed files
def read_bed(fr, chroms, bs):
    data = {k:{'+':[0]*(v//bs[k]+1), '-':[0]*(v//bs[k]+1)} for k, v in chroms.items()}
    for l in fr:
        ws = l.rstrip('\n').split('\t')
        data[ws[0]][ws[5]][int(ws[2])//bs[ws[0]]] += 1
    return data


# draw distributions
def draw(rNMPs, chroms, bs, out):
    sns.set(style='ticks')
    # get general parameters
    xmax = max(chroms.values())
    ymax = max([max(x['+'] + x['-']) for k, x in rNMPs.items() if k != 'chrM'])
    fig, axs = plt.subplots(len(chroms), figsize=(10,10))
    for i, chrom in enumerate(chroms):
        ax = axs[i]
        data = rNMPs[chrom]
        ax.scatter([x * bs[chrom] for x in range(len(data['+'])) if data['+'][x]], \
            [data['+'][x] + 10 for x in range(len(data['+'])) if data['+'][x]],\
            color='r', s=0.1)
        ax.scatter([x * bs[chrom] for x in range(len(data['-'])) if data['-'][x]], \
            [-data['-'][x] - 10 for x in range(len(data['-'])) if data['-'][x]],\
            color='b', s=0.1)
        ax.plot([chroms[chrom], chroms[chrom]], [-ymax*2, ymax*2], color='k')
        # ax.plot([0, chroms[chrom]], [0,0], color='k')
        # separate for chrM
        if chrom != 'chrM':
            ax.set_xlim((0, xmax))
            ax.set_ylim((-ymax//2, ymax//2))
        # other settings
        ax.get_xaxis().set_visible(False)
        ax.set_yticks([])
        for loc in ['top', 'bottom', 'right']:
            ax.spines[loc].set_visible(False)
        ax.set_ylabel(chrom, rotation='horizontal', labelpad=20)

        
    fig.savefig(out)


def main():
    parser = argparse.ArgumentParser(description='generate distribution of rNMP bed file in reference genome')
    parser.add_argument('bed', type=argparse.FileType('r'), help='rNMP incorporation bed file')
    parser.add_argument('fai', type=argparse.FileType('r'), help='Fasta index file for the reference genome')
    parser.add_argument('-o', default='rNMP_distribution.png', help='Output filename')
    parser.add_argument('-b', type=int, default=100000, help='Bin size (100,000)')
    parser.add_argument('-m', type=int, default=300, help='Bin size for chrM (300)')
    args = parser.parse_args()

    # generate binsizes
    binsizes = generate_bs(args.b, args.m)
    
    # read chromosome lengths
    chroms = read_fai(args.fai)

    # assign rNMPs into each bin
    rNMPs = read_bed(args.bed, chroms, binsizes)

    # draw figures
    draw(rNMPs, chroms, binsizes, args.o)

    print('Done!')

if __name__ == '__main__':
    main()
