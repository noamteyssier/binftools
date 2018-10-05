#!/usr/bin/env python3

import sys
import argparse as ap

def main():
    p = ap.ArgumentParser()
    p.add_argument('-n', '--binsize', help='length of bins')
    args = p.parse_args()

    n = int(args.binsize)
    i = 0
    current = []
    currentChrom = None
    currentPos = None


    # for i in sys.stdin:
    while True:
        try:
            chrom, pos, depth = next(sys.stdin).strip('\n').split('\t')

            if (chrom != currentChrom) and (i > 0):
                i = 0
                print(
                    '\t'.join([str(currentChrom), str(currentPos), str(int(sum(current) / n))])
                )
                currentChrom = chrom

            if i == 0:
                currentChrom = chrom

            if i%n == 0:
                print(
                    '\t'.join([str(chrom), str(pos), str(int(sum(current) / n))])
                )
                current = []

            current.append(int(depth))
            i += 1
            currentPos = pos

        except StopIteration:
            print(
                '\t'.join([str(currentChrom), str(currentPos), str(int(sum(current) / n))])
            )
            break

if __name__ == '__main__':
    main()
