#!/usr/bin/env python3

import argparse
import sys
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='Sum up match counts data')
    parser.add_argument('counts', nargs='+', type=argparse.FileType('r'), help='Matched counts data')
    parser.add_argument('-o', default=sys.stdout, type=argparse.FileType('w'), help='Output basename')
    args = parser.parse_args()

    # read data from numpy
    data = []
    df = pd.read_csv(args.counts[0], sep='\t', index_col=0)
    for f in args.counts[1:]:
        df = df.add(pd.read_csv(f, sep='\t', index_col=0))
    df.to_csv(args.o, sep='\t')

    print('Done!')

if __name__ == '__main__':
    main()
