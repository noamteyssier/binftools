#!/usr/bin/env python3

from random import choice, seed, uniform, randint
from argparse import ArgumentParser


def generate_sequences(alphabet, motifs, num_seqs, len_seqs, freq_motif):
    """generate sequences"""
    for i in range(num_seqs):

        motif = None

        # motif in sequence
        if uniform(0,1) <= freq_motif:
            motif = choice(motifs)
            len_motif = len(motif)


            pos_motif = randint(0, len_seqs - len_motif)
            seq_left = ''.join([choice(alphabet) for _ in range(pos_motif)])
            seq_right = ''.join([choice(alphabet) for _ in range(len_seqs - pos_motif - len_motif)])

            seq = seq_left + motif + seq_right

        # motif not in sequence
        else:
            seq = ''.join([choice(alphabet) for _ in range(len_seqs)])

        fasta_print(seq, motif, i)
def generate_motif(alphabet, width):
    """generate a motif of a given width"""
    return ''.join([choice(alphabet) for _ in range(width)])
def fasta_print(seq, motif, num):
    """fasta format printout"""
    print(
        '>seq_{0}_{1}\n{2}'.format(num,motif,seq)
    )
def get_args():
    """function to handle arguments; returns args object"""

    p = ArgumentParser()
    p.add_argument('-n', '--num_seqs', help='number of sequences to generate', default = 1e3)
    p.add_argument('-m', '--num_motifs', help='number of motifs to include', default = 1e2)
    p.add_argument('-l', '--len_seqs', help='length of sequences to generate', default = 1.5e2)
    p.add_argument('-f', '--freq_motif', help='frequency of finding a motif in a sequence (between 0 and 1)', default = 0.5)
    p.add_argument('-a', '--alphabet', help='[dna][rna][protein]', default = 'dna')
    p.add_argument('-w', '--motif_width', help='width of motif to add', default = 7)
    p.add_argument('-s', '--random_seed', help='optional seed for random')

    args = p.parse_args()

    assert args.alphabet in ['protein', 'dna', 'rna'], "alphabet must be protein, dna, or rna"
    assert int(args.motif_width) > 0, "motif width must be greater than 0"
    assert int(args.num_seqs) > 0, "number of sequences must be greater than 0"
    assert int(args.len_seqs) > 0, "length of sequences must be greater than 0"
    assert int(args.motif_width) <= int(args.len_seqs), "motif width must not be greater than length of sequence"

    return args

def main():
    alphabets = {
        'protein' : ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'],
        'dna' : ['A', 'C', 'T', 'G'],
        'rna' : ['A', 'C', 'U', 'G']
    }

    args = get_args()

    if args.random_seed:
        seed(int(random_seed))

    # generate motifs
    generated_motifs = [
        generate_motif(
            alphabets[args.alphabet], int(args.motif_width))
                for _ in range(int(args.num_motifs))
    ]

    # generate sequences
    generate_sequences(
        alphabets[args.alphabet],
        generated_motifs,
        num_seqs = int(args.num_seqs),
        len_seqs = int(args.len_seqs),
        freq_motif = float(args.freq_motif)
    )









if __name__ == '__main__':
    # seed(40)
    main()
