#!/usr/bin/env python3

import sys
import argparse as ap
import subprocess


def format_print(items, label, output_file):
    """method to handle print out"""
    if label:
        items.append(label)
    output_file.write('\t'.join([str(i) for i in items]) + '\n')
def collapse(bam_gen, bin_width, label, output_file):
    """main method for collapsing bam"""
    i = 0
    current = []
    currentChrom = None
    currentPos = None


    # iterate through all positions
    while True:
        try:
            chrom, pos, depth = next(bam_gen)

            # case when end of chromosome is reached before mod is zero
            if (chrom != currentChrom) and (i > 0):
                i = 0
                items = [currentChrom, currentPos, int(sum(current) / bin_width)]
                format_print(
                    items, label, output_file
                )
                currentChrom = chrom

            # reset currentChrom at new chrom
            if i == 0:
                currentChrom = chrom

            # print out mean at binsize
            if i % bin_width == 0:
                items = [chrom, pos, int(sum(current) / bin_width)]
                format_print(
                    items, label, output_file
                )
                current = []

            current.append(int(depth))
            i += 1
            currentPos = pos

        # end of loop // print what is left
        except StopIteration:
            items = [currentChrom, currentPos, int(sum(current) / bin_width)]
            format_print(
                items, label, output_file
            )
            break
def collect_depth(ifn):
    """call samtools depth and create generator"""
    p = subprocess.Popen(
        ['samtools', 'depth', '-a', ifn],
        stdout=subprocess.PIPE,
        )

    for line in iter(p.stdout.readline, b''):
        yield line.decode('ascii').strip('\n').split('\t')
def get_args():
    """method to handle arguments"""
    p = ap.ArgumentParser()
    p.add_argument('-i', '--input_file', help='bam file to collapse', required=True)
    p.add_argument('-n', '--bin_width', help='length of bins', required=False, default = 200)
    p.add_argument('-l', '--label', help = 'add label to output file as 4th column', required=False)
    p.add_argument('-o', '--output_file', help='output file to write to (default = stdout)', required=False)
    args = p.parse_args()

    # convert bin to int
    args.bin_width = int(args.bin_width)

    # write to file or stdout
    if args.output_file:
        args.output_file = open(args.output_file, 'w+')
    else:
        args.output_file = sys.stdout

    return args
def main():

    # read args
    args = get_args()

    # read bam
    bam_gen = collect_depth(args.input_file)

    # collapse bam
    collapse(bam_gen, args.bin_width, args.label, args.output_file)



if __name__ == '__main__':
    main()
