#!/usr/bin/env python3

import pandas as pd
import argparse
import gzip
import sys


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


def query_sequence(args):

    counter = 0

    for idx, record in enumerate(seqReader(args.input_sequences)):
        if record[1] == args.query_sequence:
            counter += 1
    idx += 1

    print("number of occurences : {}".format(counter))
    print("total number of reads : {}".format(idx))
    print("fraction : {}".format(counter / idx))


def sequence_hist(args):

    hist = {}

    for idx, record in enumerate(seqReader(args.input_sequences)):

        seq = record[1]

        if seq not in hist:
            hist[seq] = 0

        hist[seq] += 1


    hist = pd.DataFrame([
        {'sequence' : seq, "count": hist[seq]} for seq in hist
        ])

    hist['frequency'] = hist['count'] / idx

    hist.sort_values('frequency', inplace=True)

    hist[['count', 'frequency', 'sequence']].\
        to_csv(sys.stdout, sep="\t", index=False)


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input_sequences', required=True, type=str)
    p.add_argument('-s', '--query_sequence', required=False, type=str)
    args = p.parse_args()
    return args


if __name__ == '__main__':
    args = get_args()
    if args.query_sequence:
        query_sequence(args)
    else:
        sequence_hist(args)
