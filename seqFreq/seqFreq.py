#!/usr/bin/env python3

import argparse
import gzip


def seqReader(fn):
    """
    iterate through sequences and yield as generator
    """
    def openSeq(fn):
        if 'gz' in fn:
            return gzip.open(fn, 'rt')
        else:
            return open(fn, 'r')

    def num_iter(fn):
        if 'fastq' in fn or 'fq' in fn:
            return 4
        else:
            return 2

    n = num_iter(fn)

    with openSeq(fn) as f:
        while True:
            try:
                yield [next(f).strip('\n') for _ in range(n)]
            except StopIteration:
                break


def main(args):

    counter = 0

    for idx, record in enumerate(seqReader(args.input_sequences)):
        if record[1] == args.query_sequence:
            counter += 1
    idx += 1

    print("number of occurences : {}".format(counter))
    print("total number of reads : {}".format(idx))
    print("fraction : {}".format(counter / idx))


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input_sequences', required=True, type=str)
    p.add_argument('-s', '--query_sequence', required=True, type=str)
    args = p.parse_args()
    return args


if __name__ == '__main__':
    args = get_args()
    main(args)
