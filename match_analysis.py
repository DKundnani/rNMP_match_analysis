#!/usr/bin/env python3

import argparse
import sys

# reverse compliment
def rc(c):
    rcs = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    return rcs[c]


# Read FASTQ file to get rNMP in raw reads
# rNMP should be reverse compliment of 12th base (8 UMI and 3 barcode before rNMP)
def read_fastq(fr):
    rNMPs_raw = {}
    line_num = 1
    name = None
    for l in fr:
        if not l:
            continue
        if line_num % 4 == 1:
            assert l[0] == '@', f'Invalid read name on line {line_num} in the fastq file {fr.name}!'
            name = l[1:].split(' ')[0]
        elif line_num %4 == 2:
            rNMPs_raw[name] = rc(l[11])
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
            seq += l
    genome[name] = seq
    del genome[None]
    return genome


# find rNMPs in BED file, which is 0 based. Compare to rNMP raw reads. And return the formatted string
# Eg. chrI  15  16  .  .  + should match genome[chrI][15]
def get_rNMP_aligned(l, genome, rNMPs_raw):
    ws = l.rstrip().split('\t')
    rNMP_aligned = genome[ws[0]][int(ws[1])]
    if ws[5] == '-':
        rNMP_aligned = rc(rNMP_aligned)
    name = ws[3].split('_')[0]
    rNMP_raw = rNMPs_raw[name]
    return ws[:3] + [rNMP_raw, rNMP_aligned, ws[5]]


def main():
    parser = argparse.ArgumentParser(description='Generate bed file to check if rNMP in raw reads match with rNMP in reference genome')
    parser.add_argument('BED', nargs='+', type=argparse.FileType('r'), help='BED file for incorporated rNMPs')
    parser.add_argument('Genome', type=argparse.FileType('r'), help='Reference genome file')
    parser.add_argument('FASTQ', type=argparse.FileType('r'), help='FASTQ file for raw reads')
    parser.add_argument('-o', default='match_analysis_results', help='Output basename')
    args = parser.parse_args()

    # Process raw reads and genome
    rNMPs_raw = read_fastq(args.FASTQ)
    print('Raw reads loaded!')
    genome = read_genome(args.Genome)
    print('Genome loaded!')

    # deal with each bed file
    for bed in args.BED:
        filename = bed.name.split('/')[-1].split('.')[0]
        with open(f'{args.o}_{filename}.bed', 'w') as fw:
            for l in bed:
                data = get_rNMP_aligned(l, genome, rNMPs_raw)
                fw.write('\t'.join(data) + '\n')
        print(f'{filename} processed!')

    print('Done!')


if __name__ == '__main__':
    main()
