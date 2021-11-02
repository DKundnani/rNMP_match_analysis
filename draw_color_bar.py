#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os


def main():
    parser = argparse.ArgumentParser(description='Draw color scale only')
    parser.add_argument('tsv', type=argparse.FileType('r'), help='Tsv file with incorporated rNMP match/mismatch status')
    parser.add_argument('-o', help='Output basename')
    parser.add_argument('--palette', default='vlag', help='Seaborn palette name')
    parser.add_argument('--font_scale', default=2, type=float, help='Font scale for the marker')
    args = parser.parse_args()

    if not args.o:
        args.o = args.tsv.name.split('/')[-1].split('.')[0]

    # read data and normalize
    df = pd.read_csv(args.tsv, sep='\t')
    df['Total'] = df['A'] + df['C'] + df['G'] + df['T']
    for fea in ['A', 'C', 'G', 'T']:
        df[f'{fea}norm'] = df[fea] / df['Total']
    data = df[[f'{fea}norm' for fea in ['A', 'C', 'G', 'T']]]


    # draw
    sns.set(font_scale=args.font_scale)
    fig, ax = plt.subplots(figsize=(7,6))
    plt.subplots_adjust(top=0.98, right=0.99)
    sns.heatmap(data, vmin=0, vmax=1, annot=True, cmap=args.palette)
    ax.set_xticklabels(['A','C','G','U'])
    ax.set_yticklabels(['A','C','G','U'], rotation='horizontal')

    cbar = ax.collections[0].colorbar
    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), rotation=-90)


    # savefig
    if os.path.isdir(args.o):
        sep = '/'
    else:
        sep = '_'
    fname = args.o + sep + args.tsv.name.split('/')[-1].split('.')[0] + '_color_scale.png'
    plt.savefig(fname)


    print('Done!')

if __name__ == '__main__':
    main()
