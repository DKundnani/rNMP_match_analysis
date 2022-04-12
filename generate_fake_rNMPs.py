#!/usr/bin/env python3

import argparse
import os
import sys
import numpy.random as random

def read_fai(fr):
    res = {}
    for l in fr:
        ws = l.rstrip().split('\t')
        res[ws[0]] = int(ws[1])
    return res

def main():
    parser = argparse.ArgumentParser(description='Generate fake rNMPs from a fasta index file')
    parser.add_argument('fai', type=argparse.FileType('r'), help='Input fasta index file')
    parser.add_argument('-n', type=int, default=10000, help='Number of rNMPs (10,000)')
    parser.add_argument('--nuc', action='store_true', help='Generate on chromosomes only')
    parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output to file')
    args = parser.parse_args()

    # read chromosome length
    chroms = read_fai(args.fai)
    if args.nuc:
        del chroms['chrM']
    s = sum(chroms.values())

    # generate data
    for i in sorted(random.rand(args.n)):
        i = int(i*s*2)
        if i >= s:
            st = '-'
            i -= s
        else:
            st = '+'
        for k, v in chroms.items():
            if i < v:
                args.o.write(f'{k}\t{i}\t{i+1}\t.\t.\t{st}\n')
                break
            else:
                i -= v

    # read data
    print('Done!')

if __name__ == '__main__':
    main()
