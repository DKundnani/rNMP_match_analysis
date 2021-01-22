#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os


def main():
    parser = argparse.ArgumentParser(description='Draw heatmap for rNMP match/mismatch status')
    parser.add_argument('tsv', type=argparse.FileType('r'), help='Tsv file with incorporated rNMP match/mismatch status')
    parser.add_argument('-d', action='store_true', help='Draw heatmap for dinucleotide')
    parser.add_argument('--hide_cbar', action='store_false', help='Hide heatmap colorbar')
    parser.add_argument('-o', help='Output basename')
    args = parser.parse_args()

    if not args.o:
        args.o = args.tsv.name.split('/')[-1].split('.')[0]

    # define bases
    base = ['A','C','G','T']
    if args.d:
        temp = []
        for i in base:
            for j in base:
                temp.append(j+i)
        base = temp
    order = {x:i for i, x in enumerate(base)}
    # read data and normalize
    df = pd.read_csv(args.tsv, sep='\t')
    df = df.set_index(df.columns[0])
    df = df.transpose().reset_index()
    df = df[[df.columns[0]] + base].sort_values(by=df.columns[0], key=lambda col: col.map(order))
    df['Total'] = df.sum(axis=1)
    for fea in base:
        df[f'{fea}norm'] = df[fea] / df['Total']
    data = df[[f'{fea}norm' for fea in base]]


    # draw
    fscale = 0.9 if args.d else 2.5
    sns.set(font_scale=fscale)
    fig, ax = plt.subplots(figsize=(6,6))
    if not args.d:
        plt.subplots_adjust(top=0.99, right=0.99, left=0.11, bottom=0.11)
    else:
        plt.subplots_adjust(top=0.99, right=0.99, left=0.09, bottom=0.09)
    sns.heatmap(data, vmin=0, vmax=1, annot=True, cbar=args.hide_cbar, cmap='vlag', fmt='.0%')
    # set labels
    labelsx = ['A', 'C', 'G', 'U']
    labelsy = ['A', 'C', 'G', 'T']
    if args.d:
        tempx = []
        tempy = []
        for i in labelsx:
            for j in ['A','C','G','T']:
                tempx.append(j+i)
        for i in labelsy:
            for j in ['A','C','G','T']:
                tempy.append(j+i)
        labelsx = tempx
        labelsy = tempy
    label_size = 18 if args.d else 45
    ax.set_xticklabels(labelsx,fontsize=label_size)
    ax.set_yticklabels(labelsy, rotation='horizontal',fontsize=label_size)

    # savefig
    if os.path.isdir(args.o):
        sep = '/'
    else:
        sep = '_'
    fname = args.o + sep + args.tsv.name.split('/')[-1].split('.')[0] + '_match_analysis.png'
    plt.savefig(fname)


    print('Done!')

if __name__ == '__main__':
    main()
