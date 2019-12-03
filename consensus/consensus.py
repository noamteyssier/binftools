#!/usr/bin/env python3

import numpy as np
import argparse, sys
import pysam as ps


def parse_bed(bed_file):
    f = open(bed_file, 'r')
    while True:
        try:
            yield next(f).strip('\n').split('\t')
        except StopIteration:
            break
def RegionConsensus(aln, chrom, start, end):
    """
    get sequence consensus across a given region
    """

    start = int(start)
    end = int(end)

    region_consensus = []
    last_base=start-1
    for column in aln.pileup(chrom, start, end, truncate=True):
        if column.pos != last_base+1:
            [region_consensus.append('N') for _ in range(column.pos - start)]

        try:
            seq_column = np.array([b.upper() for b in column.get_query_sequences()])

            bases, counts = np.unique(seq_column, return_counts=True)

            consensus = bases[np.argmax(counts)]

            region_consensus.append(consensus)

            last_base = column.pos

        except AssertionError:
            break

    # case where full region is missing
    if len(region_consensus) == 0:
        size = end - start
        region_consensus = ['N'] * size

    return ''.join(region_consensus)
def WriteFasta(consensus, label, file_out):
    """
    write fasta given consensus sequence, label, and output file
    """
    file_out.write(
        '>{}\n{}\n'.format(label, consensus)
    )

def main(args):

    output_file = open(args.output_file, 'w+') if args.output_file else sys.stdout
    bed_region = next(parse_bed(args.bed_file))
    aln = ps.AlignmentFile(args.input_file)

    for i, bed_region in enumerate(parse_bed(args.bed_file)):
        chrom, start, end = bed_region[:3]
        seq = RegionConsensus(aln, chrom, start, end)


        try:
            WriteFasta(seq, bed_region[3], output_file)
        except IndexError:
            WriteFasta(seq, i, output_file)


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("-i", '--input_file',
        help='input bam to generate consensus sequence over regions'
        )
    p.add_argument("-b", '--bed_file',
        help='bed file of regions'
        )
    p.add_argument('-r', '--region',
         help='standalone region to generate consensus'
         )
    p.add_argument('-o', '--output_file',
         help='where to write output fasta file (default = stdout)'
         )
    p.add_argument('-l', '--label_column',
        help='which column of bed file to use as label (default = 4 // if none found will label numerically)'
        )
    args = p.parse_args()

    if (not args.bed_file) and (not args.region):
        sys.exit('\nError : requires either bed file or region\n')

    return args
if __name__ == '__main__':
    args = get_args()
    main(args)
