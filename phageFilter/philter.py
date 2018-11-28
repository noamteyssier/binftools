#!/usr/bin/env python3

import argparse
import gzip
import re
import subprocess
import sys

def write_fastq(fq, data):
    """write lines of fastq to file"""
    [fq.write(line.decode('ascii')) for line in data]
def cut_fastq(fq, pattern):
    """trim length of pattern from read"""
    fq[1] = fq[1].decode('ascii')[len(pattern):].encode()
    fq[3] = fq[3].decode('ascii')[len(pattern):].encode()
    return fq
def main():
    p = argparse.ArgumentParser()
    p.add_argument('-1', '--forward_fq', help='filename of forward fastq', required=True)
    p.add_argument('-2', '--reverse_fq', help='filename of reverse fastq', required=True)
    p.add_argument('-o', '--output_base', help='base name of output (philtered_$base)', required=True)
    p.add_argument('-z', '--gzip', action='store_true', help='compress file after completion', required=False)
    args = p.parse_args()

    f_pattern = 'CGAATTCAGTGGTTGGTGCTGTAGGAGCA'
    r_pattern = 'AAGCTTGAGGCCATGGCATATGC'

    forward_fq = gzip.open(args.forward_fq, 'r')
    reverse_fq = gzip.open(args.reverse_fq, 'r')

    reg_forward = re.compile('^' + f_pattern)
    reg_reverse = re.compile('^' + r_pattern)


    f_out = open('philtered_{0}_R1.fq'.format(args.output_base), 'w+')
    r_out = open('philtered_{0}_R2.fq'.format(args.output_base), 'w+')


    pairs_kept = 0
    pairs_total = 0

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

    # zip both files
    if args.gzip:
        subprocess.Popen(
            "gzip philtered_{0}_R1.fq".format(args.output_base),
            shell=True
        ).wait()
        subprocess.Popen(
            "gzip philtered_{0}_R2.fq".format(args.output_base),
            shell=True
        ).wait()

    print('Pairs Kept   : {0}'.format(pairs_kept))
    print('Total Pairs  : {0}'.format(pairs_total))
    print('Percent Kept : {0}'.format(pairs_kept/pairs_total))


if __name__ == '__main__':
    main()
