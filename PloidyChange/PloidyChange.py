#!/usr/bin/env python3

import argparse, sys
import pysam as ps

def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("-i", '--input_vcf',
        default=sys.stdin,
        help='input_file (default = stdin)'
        )
    p.add_argument("-o", '--output_vcf',
        help='input_file (default = stdout)'
        )
    p.add_argument("-p", '--ploidy',
        default=2, type=int,
        help='ploidy to convert to (default=2)'
        )
    args = p.parse_args()
    return args
def write_header(vcf, fileout):
    print(
        vcf.header.__str__().strip('\n'),
        file=fileout
        )
def write_record(record, fileout):
    print(
        record,
        file=fileout
    )
def split_record(record):
    return record.__str__().strip('\n').split('\t')
def convert_gt(gt_info, ploidy=2):
    return ['|'.join([gt] * ploidy) for gt in gt_info]
def convert_record(record, ploidy=2):
    info = split_record(record)
    meta = info[:9]
    fixed_gt = convert_gt(info[9:], ploidy)
    return '\t'.join(meta + fixed_gt)
def main(args):
    vcf = ps.VariantFile(args.input_vcf)
    fileout = sys.stdout if not args.output_vcf else open(args.output_vcf, 'w+')

    write_header(vcf, fileout)
    for record in vcf:
        write_record(
            convert_record(record, args.ploidy),
            fileout
            )

if __name__ == '__main__':
    args = get_args()
    main(args)
