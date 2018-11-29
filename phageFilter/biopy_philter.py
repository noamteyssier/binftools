#!/usr/bin/env python3

import argparse
import re
import sys
import time
import os
from Bio import SeqIO
from Bio.Seq import Seq
from gzip import open as gzopen


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
def status_update(pairs_kept, pairs_total, start_time):
    eprint('{0:10.4f}s:{1:10.4f}% kept of {2} pairs processed'.format(time.time() - start_time, pairs_kept/pairs_total, pairs_total))
def final_printout(pairs_kept, pairs_total, start_time):
    eprint('-------')
    eprint('Pairs Kept   : {0}'.format(pairs_kept))
    eprint('Total Pairs  : {0}'.format(pairs_total))
    eprint('Percent Kept : {0}'.format(pairs_kept/pairs_total))
    eprint('Time Elapsed : {0}'.format(time.time() - start_time))
    eprint('-------\n')
def cut_fastq(fq, pattern):
    """trim length of pattern from read"""

    l = len(pattern)

    new_quality = fq.letter_annotations['phred_quality'][l:]

    fq.letter_annotations = {}
    fq.seq = Seq(
        str(fq.seq)[l:]
    )
    fq.letter_annotations = {
        'phred_quality' : new_quality
    }

    return fq

def write_fastq(fn, fq):
    SeqIO.write(
        fq, fn, 'fastq'
    )



def main():
    p = argparse.ArgumentParser()
    p.add_argument('-1', '--forward_fq', help='filename of forward fastq', required=True)
    p.add_argument('-2', '--reverse_fq', help='filename of reverse fastq', required=True)
    p.add_argument('-o', '--output_base', help='base name of output (philtered_$base)', required=True)
    p.add_argument('-z', '--gzip', action='store_true', help='compress file after completion', required=False)
    args = p.parse_args()

    start_time = time.time()

    # assign output names
    f_ofn = 'philtered_{0}_R1.fq.gz'.format(args.output_base)
    r_ofn = 'philtered_{0}_R2.fq.gz'.format(args.output_base)

    # patterns to trim
    f_pattern = 'CGAATTCAGTGGTTGGTGCTGTAGGAGCA'
    r_pattern = 'AAGCTTGAGGCCATGGCATATGC'

    # input files
    forward_fq = SeqIO.parse(gzopen(args.forward_fq, 'rt'), "fastq")
    reverse_fq = SeqIO.parse(gzopen(args.reverse_fq, 'rt'), "fastq")

    # output files
    f_out = gzopen(f_ofn, 'wt')
    r_out = gzopen(r_ofn, 'wt')

    # regex of patterns
    reg_forward = re.compile('^' + f_pattern)
    reg_reverse = re.compile('^' + r_pattern)

    # output totals
    pairs_kept = 0
    pairs_total = 0


    # with bgzf.BgzfWriter("test.fastq.bgz", "wb") as outgz:
    #     SeqIO.write(sequences=records, handle=outgz, format="fastq")

    while True:
        pairs_total += 1
        try:
            f,r = next(forward_fq), next(reverse_fq)

            if pairs_total % 100000 == 0:
                status_update(pairs_kept, pairs_total, start_time)


            if reg_forward.match(str(f.seq)) and reg_reverse.match(str(r.seq)):
                cut_fastq(f, f_pattern)
                cut_fastq(r, r_pattern)
                write_fastq(f_out, f)
                write_fastq(r_out, r)
                pairs_kept += 1


        except StopIteration:
            [i.close() for i in [f_out, r_out]]
            final_printout(pairs_kept, pairs_total, start_time)
            break



    # # loop over all reads
    # while True:
    #     pairs_total += 1
    #     try:
    #         f = [next(forward_fq) for _ in range(4)]
    #         r = [next(reverse_fq) for _ in range(4)]
    #
    #         if pairs_total % 100000 == 0:
    #             status_update(pairs_kept, pairs_total, start_time)
    #
    #         # only proceed if both pairs match regex
    #         if reg_forward.match(f[1].decode('ascii')) and reg_reverse.match(r[1].decode('ascii')):
    #             f = cut_fastq(f, f_pattern)
    #             r = cut_fastq(r, r_pattern)
    #             write_fastq(f_out, f)
    #             write_fastq(r_out, r)
    #             pairs_kept += 1

    #     except StopIteration:
    #         f_out.close()
    #         r_out.close()
    #         if args.gzip:
    #             gzip_file(f_ofn)
    #             gzip_file(r_ofn)
    #         break
    #     except OSError:
    #         subprocess.Popen(
    #             "rm philtered_{0}_R1.fq".format(args.output_base),
    #             shell=True
    #         )
    #         subprocess.Popen(
    #             "rm philtered_{0}_R2.fq".format(args.output_base),
    #             shell=True
    #         )
    #         sys.exit('ERROR : input files are not gzip')
    #
    # final_printout(pairs_kept, pairs_total, start_time)
    #

if __name__ == '__main__':
    main()
