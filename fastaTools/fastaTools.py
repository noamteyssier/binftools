#!/usr/bin/env python3

import sys, argparse
import numpy as np

def read_fasta(fn):
    if fn != sys.stdin:
        f = open(fn, 'rt')
    else:
        f = fn

    while True:
        try:
            header, seq = [next(f) for _ in range(2)]
            yield header, seq
        except StopIteration:
            f.close()
            break
def write_fasta(header, seq):
    sys.stdout.write("{}{}".format(header, seq))
def get_args():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input_fn', help='fasta file to process', default=sys.stdin)
    p.add_argument('-s', '--subsample', help='percentage to subsample', default=None, type=float)
    args = p.parse_args()
    return args

def subsample(args):
    for header, seq in read_fasta(args.input_fn):
        if np.random.random() <= args.subsample:
            write_fasta(header, seq)
def main():
    args = get_args()
    if args.subsample:
        subsample(args)


    pass

if __name__ == '__main__':
    main()
