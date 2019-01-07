#!/usr/bin/env python3

import numpy as np
from random import choice, sample, seed, uniform, randint
from argparse import ArgumentParser
from itertools import product, chain

def generate_snps(alphabet, numSNP):
    """return a list of snps in the variant format [major, minor]"""
    return [sample(alphabet,2) for i in range(numSNP)]

def generate_snp_positions(numSNP, lenSequences):
    """return a list of positions to create snps"""
    return sorted(sample(range(lenSequences), numSNP))

def generate_haplotypes(alphabet, lenSequences, numSNP, numHaplotypes):
    """create the base amplicon to create haplotypes within"""
    snps = generate_snps(alphabet, numSNP)
    snp_positions = generate_snp_positions(numSNP, lenSequences)
    snp_combos = list(product(*snps))
    base_amplicon = [choice(alphabet) for _ in range(lenSequences)]
    haplotypes_to_use = sample(range(len(snp_combos)), numHaplotypes)

    for h in range(numHaplotypes):
        current_haplotype = base_amplicon.copy()
        for i in range(len(snp_positions)):
            current_haplotype[snp_positions[i]] = snp_combos[haplotypes_to_use[h]][i]
        yield current_haplotype

def introduce_error(alphabet, seq, errorRate):
    """introduce a random base relative to the error rate"""
    new_seq = seq.copy()
    for i in range(len(seq)):
        if uniform(0,1) <= errorRate:
            new_seq[i] = choice(alphabet)
    return new_seq

def generate_sequences(alphabet, numSequences, numHaplotypes, numSNP, errorRate, lenSequences):
    """main loop : calls haplotypes, introduces errors, and sends to fasta printout"""
    possible_haplotypes = list(generate_haplotypes(alphabet, lenSequences, numSNP, numHaplotypes))
    for sequence_number in range(numSequences):
        base_haplotype = choice(possible_haplotypes).copy()
        haplotype_number = possible_haplotypes.index(base_haplotype)
        haplotype = introduce_error(alphabet, base_haplotype, errorRate)
        fasta_print(sequence_number, haplotype_number, haplotype)


def fasta_print(sequence_number, haplotype_number, haplotype):
    """print in fasta format with relevant information"""
    print(
        '>seq{0}_hap{1}\n{2}'.format(sequence_number, haplotype_number, ''.join(haplotype))
    )


def get_args():
    """function to handle arguments; returns args object"""

    p = ArgumentParser()
    p.add_argument('-n', '--num_seqs', help='number of sequences to generate', default = 1e3)
    p.add_argument('-k', '--num_haplotypes', help='number of haplotypes to include', default = 2)
    p.add_argument('-m', '--num_snps', help='number of snp positions to use', default = 10)
    p.add_argument('-l', '--len_seqs', help='length of sequences to generate', default = 1.5e2)
    p.add_argument('-e', '--error_rate', help='frequency of introducing error into sequence (between 0 and 1)', default = 1e-3)
    p.add_argument('-a', '--alphabet', help='[dna][rna][protein]', default = 'dna')
    p.add_argument('-s', '--random_seed', help='optional seed for random')

    args = p.parse_args()

    assert args.alphabet in ['protein', 'dna', 'rna'], "alphabet must be protein, dna, or rna"
    assert int(args.num_haplotypes) > 0, "number of haplotypes must be greater than 0"
    assert int(args.num_seqs) > 0, "number of sequences must be greater than 0"
    assert int(args.len_seqs) > 0, "length of sequences must be greater than 0"
    assert int(args.num_snps) >= int(args.num_haplotypes), "number of snps must be greater than or equal to number of haplotypes"

    return args

def main():
    alphabets = {
        'protein' : ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'],
        'dna' : ['A', 'C', 'T', 'G'],
        'rna' : ['A', 'C', 'U', 'G']
    }

    args = get_args()

    alphabet = alphabets[args.alphabet]

    numSequences = int(args.num_seqs)
    numHaplotypes = int(args.num_haplotypes)
    numSNP = int(args.num_snps)
    errorRate = float(args.error_rate)
    lenSequences = int(args.len_seqs)

    if args.random_seed:
        seed(args.random_seed)


    generate_sequences(
        alphabet,
        numSequences,
        numHaplotypes,
        numSNP,
        errorRate,
        lenSequences
    )


if __name__ == '__main__':
    main()
