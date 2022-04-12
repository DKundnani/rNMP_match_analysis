#!/usr/bin/env python3

import argparse
import sys

# reverse compliment
def rc(s):
    rcs = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    out = ''
    for c in s:
        out = rcs[c] + out
    return out


# Read FASTQ file to get rNMP in raw reads
# rNMP should be reverse compliment of 12th base (8 UMI and 3 barcode before rNMP)
def read_fastq(fr, deviation, n):
    rNMPs_raw = {}
    line_num = 1
    name = None
    for l in fr:
        l = l.rstrip('\n')
        if not l:
            continue
        if line_num % 4 == 1:
            assert l[0] == '@', f'Invalid read name on line {line_num} in the fastq file {fr.name}!'
            name = l[1:].split(' ')[0]
        elif line_num % 4 == 2:
            rNMPs_raw[name] = rc(l[deviation:min(len(l), deviation+n)])
        line_num += 1
    return rNMPs_raw


# read genome
def read_genome(fr):
    genome = {}
    name = None
    seq = None
    for l in fr:
        l = l.rstrip()
        if not l:
            continue
        if l[0] == '>':
            genome[name] = seq
            name = l[1:]
            seq = ''
        else:
            seq += l.upper()
    genome[name] = seq
    del genome[None]
    return genome


# find rNMPs in BED file, which is 0 based. Compare to rNMP raw reads. And return the formatted string
# Eg. chrI  15  16  .  .  + should match genome[chrI][15]
def get_rNMP_aligned(l, genome, rNMPs_raw, n):
    ws = l.rstrip().split('\t')
    loc = int(ws[1])
    # get base
    if ws[5] == '+':
        rNMP_aligned = genome[ws[0]][max(0, loc+1-n):loc+1]
    else:
        rNMP_aligned = rc(genome[ws[0]][loc:min(loc+n, len(genome[ws[0]]))])
    name = ws[3].split('_')[0]
    rNMP_raw = rNMPs_raw[name]
    return [ws[0], ws[1], ws[2], rNMP_raw, rNMP_aligned, ws[5]]


def main():
    parser = argparse.ArgumentParser(description='Generate bed file to check if rNMP in raw reads' + \
                                                    'match with rNMP in reference genome')
    parser.add_argument('BED', nargs='+', type=argparse.FileType('r'), help='BED file for incorporated rNMPs')
    parser.add_argument('Genome', type=argparse.FileType('r'), help='Reference genome file')
    parser.add_argument('FASTQ', type=argparse.FileType('r'), help='FASTQ file for raw reads')
    parser.add_argument('-d', type=int, default=11, help='Deviation to the read ends. '+ \
                                                        ' (11nt for ribose-seq, 0 for emriboseq), (11)')
    parser.add_argument('-n', type=int, default=1, help='Number of nucleotides to compare, the ' + \
                                                        'last one is rNMP. n=1 means rNMP only (1)')
    parser.add_argument('-o', default='match_analysis_results', help='Output basename')
    args = parser.parse_args()

    # Process raw reads and genome
    rNMPs_raw = read_fastq(args.FASTQ, args.d, args.n)
    print('Raw reads loaded!')
    genome = read_genome(args.Genome)
    print('Genome loaded!')

    # deal with each bed file
    for bed in args.BED:
        filename = bed.name.split('/')[-1].split('.')[0]
        with open(f'{args.o}_{filename}.bed', 'w') as fw:
            for l in bed:
                data = get_rNMP_aligned(l, genome, rNMPs_raw, args.n)
                if data:
                    fw.write('\t'.join(data) + '\n')
        print(f'{filename} processed!')

    print('Done!')


if __name__ == '__main__':
    main()
