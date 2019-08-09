#!/usr/bin/env python3

import sys, argparse, gzip


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("-i",
        "--input_positions",
        help='file containing at least two tabs, [chrom, position, etc...]',
        required=True
        )
    p.add_argument("-n",
        "--window_size",
        help='expansion side on either side of position (n=50 means 100bp window)',
        required=False,
        default=50,
        type=int
        )
    p.add_argument('--no_label',
        action='store_false',
        help='Do not label expanded bed'
        )
    args = p.parse_args()
    return args
def iter_positions(filename):
    """iterate through positions"""
    f = open(filename, 'r') if '.gz' not in filename else gzip.open(filename, 'rt')
    while True:
        try:
            line = next(f)
            if '#' not in line:
                yield line
        except StopIteration:
            break
def add_window(position, size):
    """expand position"""
    position = int(position)
    start = position - size
    stop = position + size
    if start < 0 :
        start = 0
    return start, stop
def write_out(chrom, start, stop, count):
    """write bed to stdout"""
    sys.stdout.write('\t'.join(
        [str(i) for i in [chrom, start, stop, count]]
    ) + '\n')


def main(args):

    count = 0
    for line in iter_positions(args.input_positions):
        vars = line.split('\t')
        start, stop = add_window(vars[1], args.window_size)
        write_out(vars[0], start, stop, count)
        count += 1


if __name__ == '__main__':
    args = get_args()
    main(args)
