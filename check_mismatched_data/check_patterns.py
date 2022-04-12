#!/usr/bin/env python3

import argparse
import sys
import itertools as it
import numpy as np


# read bed file
def read_fa(fr, k):
    res = {}
    for l in fr:
        if l[0] == '>':
            continue
        l = l.rstrip('\n').upper()
        for i in range(len(l)+1-k):
            kmer = l[i:i+k]
            if kmer not in res:
                res[kmer] = 0
            res[kmer] += 1
    freqs = [(k,v/sum(res.values())) for k,v in res.items()]
    return sorted(freqs, key=lambda x: -x[1])


# output to files
def output(freqs, fw):
    n = max(len(x) for x in freqs.values())
    fw.write(f'Order\t' + '\t'.join([f'{x}_pattern\t{x}_freq' for x in freqs.keys()]) + '\n')
    for i in range(n):
        fw.write(f'{i+1}')
        for k, v in freqs.items():
            if len(v) < i+1:
                fw.write('\t\t')
            else:
                fw.write(f'\t{v[i][0]}\t{v[i][1]:.4f}')
        fw.write('\n')


def main():
    parser = argparse.ArgumentParser(description='Check composition at each location')
    parser.add_argument('fas', type=argparse.FileType('r'), nargs='+', help='Fasta files')
    parser.add_argument('-k', type=int, default=1, help='K-mer length, default=1')
    parser.add_argument('-o', default=sys.stdout, type=argparse.FileType('w'), help='Output to file')
    args = parser.parse_args()

    # get kmer frequencies for each fasta files
    freqs = {fa.name.split('/')[-1].split('.')[0]:read_fa(fa, args.k) for fa in args.fas}

    # output to file
    output(freqs, args.o)

    print('Done!')

if __name__ == '__main__':
    main()
