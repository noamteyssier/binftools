#!/usr/bin/env python3

import sys
import argparse as ap

def main():
    p = ap.ArgumentParser()
    p.add_argument('-n', '--binsize', help='length of bins', required=True)
    args = p.parse_args()

    n = int(args.binsize)
    i = 0
    current = []
    currentChrom = None
    currentPos = None


    # iterate through all positions
    while True:
        try:
            chrom, pos, depth = next(sys.stdin).strip('\n').split('\t')

            # case when end of chromosome is reached before mod is zero
            if (chrom != currentChrom) and (i > 0):
                i = 0
                print(
                    '\t'.join([str(currentChrom), str(currentPos), str(int(sum(current) / n))])
                )
                currentChrom = chrom

            # reset currentChrom at new chrom
            if i == 0:
                currentChrom = chrom

            # print out mean at binsize
            if i%n == 0:
                print(
                    '\t'.join([str(chrom), str(pos), str(int(sum(current) / n))])
                )
                current = []

            current.append(int(depth))
            i += 1
            currentPos = pos

        # end of loop // print what is left
        except StopIteration:
            print(
                '\t'.join([str(currentChrom), str(currentPos), str(int(sum(current) / n))])
            )
            break

if __name__ == '__main__':
    main()
