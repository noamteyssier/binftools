#!/usr/bin/env python3

import numpy as np
import pandas as pd
from tqdm import *
import argparse, sys, gzip


class PileupToAlleleFrequency:
    """
    Convert a pileup to an allele frequency matrix
    """
    def __init__(self, pileup_file, minimum_depth=5):
        self.fn = pileup_file
        self.minimum_depth = minimum_depth

        self.alphabet = ['A', 'C', 'T', 'G']
        self.loci = []
        self.frequencies = []

        self.CountMatrix = np.array([])
        self.LociSums = np.array([])
        self.FrequencyMatrix = np.array([])
        self.FrequencyTable = pd.DataFrame()
    def IterPile(self):
        """
        generator for pileup
        """
        f = open(self.fn, 'r') if '.gz' not in self.fn else gzip.open(self.fn, 'rt')
        while True:
            try:
                yield next(f).strip('\n').split('\t')
            except StopIteration:
                break
    def PileToArray(self, pile):
        """
        Keep only ACTG, return as np.array
        """
        return np.array([
            self.alphabet.index(c) for c in pile.upper() if c in self.alphabet
                ])
    def PileCounts(self, pile):
        """
        Number of each unique allele in pileup
        """
        alleles, counts = np.unique(
            self.PileToArray(pile),
            return_counts=True
            )
        return alleles, counts
    def AppendLocusCounts(self, chrom, pos, alleles, counts):
        """
        Index loci
        Increment allele count
        """
        locus = ':'.join([chrom, pos])

        if locus not in self.loci:
            self.loci.append(locus)
            self.frequencies.append(np.zeros(4))


        for i, base in enumerate(alleles):
            self.frequencies[self.loci.index(locus)][base] += counts[i]
    def BuildCountMatrix(self):
        """
        Iterate over pileups and append to locus count matrix
        """
        for row in tqdm(self.IterPile(), desc='iterating pileup'):
            chrom = row[0]
            pos = row[1]
            depth = row[3]
            pile = row[4]

            if int(depth) < self.minimum_depth:
                continue

            # get unique alleles and respective counts
            alleles, counts = self.PileCounts(pile)

            # add to locus counts
            self.AppendLocusCounts(chrom, pos, alleles, counts)

        self.CountMatrix = np.vstack(self.frequencies)
    def BuildFrequencyMatrix(self):
        """
        Convert count matrix to frequency matrix
        """
        self.LociSums = self.CountMatrix.sum(axis=1)
        self.FrequencyMatrix = (
            self.CountMatrix.T / self.LociSums
            ).T
    def fit(self):
        """
        Build count vector for all loci
        convert loci count matrix to frequencies
        return as pandas dataframe
        """
        self.BuildCountMatrix()
        self.BuildFrequencyMatrix()

        self.FrequencyTable = pd.DataFrame(
            self.FrequencyMatrix,
            columns = self.alphabet
            )
        self.FrequencyTable['locus'] = self.loci
        self.FrequencyTable['depth'] = self.LociSums.astype(int)

        self.FrequencyTable.sort_values('locus', inplace=True)

        return self.FrequencyTable

def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("-i", '--pileup_fn',
        required=True,
        type=str,
        help='input_file'
        )
    p.add_argument("-d", '--minimum_depth',
        required=False,
        default=5,
        type=int,
        help='minimum depth to consider a pileup'
        )

    args = p.parse_args()
    return args
def main(args):
    p2af = PileupToAlleleFrequency(args.pileup_fn, minimum_depth=args.minimum_depth)
    fm = p2af.fit()
    fm.to_csv(sys.stdout, sep="\t", index=False)


if __name__ == '__main__':
    args = get_args()
    main(args)
