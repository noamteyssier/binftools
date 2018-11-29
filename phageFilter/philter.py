#!/usr/bin/env python3

import argparse
import gzip
import re
import subprocess
import sys
import time
import shutil
import os

def write_fastq(fq, data):
    """write lines of fastq to file"""
    [fq.write(line.decode('ascii')) for line in data]


def cut_fastq(fq, pattern):
    """trim length of pattern from read"""
    fq[1] = fq[1].decode('ascii')[len(pattern):].encode()   # sequence trim
    fq[3] = fq[3].decode('ascii')[len(pattern):].encode()   # quality trim
    return fq

def gzip_file(fn):
    """gzip an existing file"""
    with open(fn, 'rb') as f_in:
        with gzip.open(fn+'.gz', 'wb+') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(fn)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('-1', '--forward_fq', help='filename of forward fastq', required=True)
    p.add_argument('-2', '--reverse_fq', help='filename of reverse fastq', required=True)
    p.add_argument('-o', '--output_base', help='base name of output (philtered_$base)', required=True)
    p.add_argument('-z', '--gzip', action='store_true', help='compress file after completion', required=False)
    args = p.parse_args()

    start_time = time.time()

    # assign output names
    f_ofn = 'philtered_{0}_R1.fq'.format(args.output_base)
    r_ofn = 'philtered_{0}_R2.fq'.format(args.output_base)

    # patterns to trim
    f_pattern = 'CGAATTCAGTGGTTGGTGCTGTAGGAGCA'
    r_pattern = 'AAGCTTGAGGCCATGGCATATGC'

    # input files
    forward_fq = gzip.open(args.forward_fq, 'r')
    reverse_fq = gzip.open(args.reverse_fq, 'r')

    # output files
    f_out = open(f_ofn, 'w+')
    r_out = open(r_ofn, 'w+')

    # regex of patterns
    reg_forward = re.compile('^' + f_pattern)
    reg_reverse = re.compile('^' + r_pattern)

    # output totals
    pairs_kept = 0
    pairs_total = 0

    # loop over all reads
    while True:
        pairs_total += 1
        try:
            f = [next(forward_fq) for _ in range(4)]
            r = [next(reverse_fq) for _ in range(4)]

            # only proceed if both pairs match regex
            if reg_forward.match(f[1].decode('ascii')) and reg_reverse.match(r[1].decode('ascii')):
                f = cut_fastq(f, f_pattern)
                r = cut_fastq(r, r_pattern)
                write_fastq(f_out, f)
                write_fastq(r_out, r)
                pairs_kept += 1
        except StopIteration:
            f_out.close()
            r_out.close()
            if args.gzip:
                gzip_file(f_ofn)
                gzip_file(r_ofn)
            break
        except OSError:
            subprocess.Popen(
                "rm philtered_{0}_R1.fq".format(args.output_base),
                shell=True
            )
            subprocess.Popen(
                "rm philtered_{0}_R2.fq".format(args.output_base),
                shell=True
            )
            sys.exit('ERROR : input files are not gzip')

    sys.exit()


    print('-------')
    print('Pairs Kept   : {0}'.format(pairs_kept))
    print('Total Pairs  : {0}'.format(pairs_total))
    print('Percent Kept : {0}'.format(pairs_kept/pairs_total))
    print('Time Elapsed : {0}'.format(time.time() - start_time))
    print('-------\n')


if __name__ == '__main__':
    main()
