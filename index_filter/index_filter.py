#!/usr/bin/env python3

import argparse, sys, gzip

def get_index(record):
    name, index = record[0].split(' ')
    return index.split(':')[-1]
def read_header(fastq):
    indices = {}
    total_reads = 0
    f = gzip.open(fastq, 'rt')
    while True:
        try:
            record = [next(f) for _ in range(4)]
            index = get_index(record)

            if index not in indices:
                indices[index] = 0
            indices[index] += 1
            total_reads += 1
        except StopIteration:
            f.close()
            max_index = max(indices, key=indices.get)
            return max_index, indices[max_index], total_reads
def write_record(f, record):
    [f.write(r) for r in record]
def select_index(fastq_in, fastq_out, max_index):
    f_in = gzip.open(fastq_in, 'rt')
    f_out = gzip.open(fastq_out, 'wt')

    while True:
        try:
            record = [next(f_in) for _ in range(4)]
            if get_index(record) == max_index:
                write_record(f_out, record)

        except StopIteration:
            f_in.close()
            f_out.close()
            return




def get_args():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input1', help='read 1 of fastq', required=True)
    p.add_argument('-I', '--input2', help='read 2 of fastq', required=True)
    p.add_argument('-o', '--output1', help='output 1 of fastq', required=True)
    p.add_argument('-O', '--output2', help='output 2 of fastq', required=True)
    args = p.parse_args()
    return args

def main():
    args = get_args()
    max_index, idx_reads, total_reads = read_header(args.input1)

    select_index(args.input1, args.output1, max_index)
    select_index(args.input2, args.output2, max_index)
    sys.stderr.write(
        '{} / {} reads kept ({:.3f}% passing filter)\n'.format(
            idx_reads, total_reads, idx_reads/total_reads
        ))

if __name__ == '__main__':
    main()
