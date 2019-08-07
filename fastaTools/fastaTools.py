#!/usr/bin/env python3

import sys, argparse
import numpy as np
import gzip

def read_fasta(fn):
    if fn != sys.stdin:
        f = open(fn, 'rt') if '.gz' != fn[-3:] else gzip.open(fn, 'rt')
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
def subsample(args):
    for header, seq in read_fasta(args.input_fn):
        if np.random.random() <= args.subsample:
            write_fasta(header, seq)
def translate(args):
    def chunker(seq, size):
        return (seq[pos:pos + size] for pos in np.arange(0, len(seq), size))
    translation = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

    for header, seq in read_fasta(args.input_fn):
        if len(seq) >= 3:
            to_translate = seq[args.offset:].strip('\n')
            aa = ''.join(
                    [translation[codon] for codon in chunker(to_translate, 3) if len(codon) == 3]
                    ) + '\n'
            write_fasta(header, aa)
def main(args):
    if args.subsample:
        subsample(args)
        return
    if args.translate:
        translate(args)


    pass
def get_args():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input_fn', help='fasta file to process', default=sys.stdin)
    p.add_argument('-s', '--subsample', help='percentage to subsample', default=None, type=float)
    p.add_argument('-t', '--translate', help='translate to amino acid', action='store_true')
    p.add_argument('--offset', help='translation offset', default=0, type=int)
    args = p.parse_args()
    return args
if __name__ == '__main__':
    args = get_args()
    main(args)
