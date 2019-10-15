#!/usr/bin/env python3

import numpy as np
from nltk import ngrams
from sklearn.decomposition import PCA
from scipy.spatial.distance import *

import seaborn as sns
import matplotlib.pyplot as plt

import gzip
import sys

class EmbedSeq(object):

    def __init__(self, **kwargs):
        self.expected_length = 235

        self.sequence_table = {}

        self.sequence_id = {}

        self.matrix = np.array([])

        self.alphabet = {
            'A': 0,
            'C': 1,
            'T': 2,
            'G': 3
            }

        self.drop_n = True

        self.drop_length = True

    def iter_seq(self, fn):

        count = 0

        for record in self.sequence_reader(fn):

            seq = record[1]

            if self.drop_length and len(seq) != self.expected_length:
                continue

            if self.drop_n and 'N' in seq:
                continue

            if seq not in self.sequence_table:
                self.sequence_table[seq] = 1
                self.sequence_id[count] = seq
                count += 1
                yield seq

            else:
                self.sequence_table[seq] += 1

    def fit(self, fn):

        self.__fit__(fn)

    def transform(self):

        return self.matrix

    def fit_transform(self, fn):

        self.fit(fn)

        return self.transform()

    def sequence_reader(self, fn):
        def open_file(fn):
            if 'gz' in fn:
                return gzip.open(fn, 'rt')
            else:
                return open(fn, 'r')
        def is_fastq(fn):
            return ('fq' in fn) or ('fastq' in fn)

        f = open_file(fn)
        steps = 4 if is_fastq(fn) else 2
        while True:
            try:
                yield [next(f).strip('\n') for _ in range(steps)]
            except StopIteration:
                break


class seq_to_onehot(EmbedSeq):
    """
    Seq to One Hot Encoding
    """
    def __init__(self, **kwargs):

        EmbedSeq.__init__(self, **kwargs)

        self.__dict__.update(kwargs)

    def __fit__(self, fn):

        self.matrix = np.stack(
            [self.convert(s) for s in self.iter_seq(fn)]
            )

    def convert(self, seq):

        arr = np.zeros(
            (self.expected_length, len(self.alphabet))
            )

        for i, base in enumerate(seq):
            arr[i][self.alphabet[base]] = 1

        return arr


class seq_to_bag(EmbedSeq):
    """
    Sequences to bag of words style ngram encoding
    """
    def __init__(self, **kwargs):

        EmbedSeq.__init__(self, **kwargs)

        self.ngram_table = {}
        self.ngram_size = 7
        self.drop_n = False
        self.seq_id = 0

        self.__dict__.update(kwargs)

    def __fit__(self, fn):

        self.matrix = np.stack(
            [i for i in self.build_bag(fn)]
            )

    def build_bag(self, fn):
        samples = [self.convert(s) for s in self.iter_seq(fn)]
        for s in samples:
            arr = np.zeros(self.seq_id)
            arr[s] = 1
            yield arr

    def convert(self, seq):
        id_arr = []

        for i in ngrams(seq, self.ngram_size):
            ng = ''.join(i)

            if ng not in self.ngram_table:
                self.ngram_table[ng] = self.seq_id
                self.seq_id += 1

            id_arr.append(
                self.ngram_table[ng]
                )

        return np.array(id_arr)


def plot_pca(mat):
    pca = PCA(n_components=10)
    pcs = pca.fit_transform(mat)
    sns.scatterplot(
        x = pcs[:,0],
        y = pcs[:,1]
        )
    plt.show()

def to_onehot(fasta):
    oh = seq_to_onehot()
    mat = oh.fit_transform(fasta)
    flattened = np.stack([i.ravel() for i in mat])
    plot_pca(flattened)

def to_ngrams(fasta):
    ng = seq_to_bag(ngram_size = 10)
    mat = ng.fit_transform(fasta)

    plot_pca(mat)

if __name__ == '__main__':
    fasta = "downsampled_ts.fa"
    to_onehot(fasta)
    to_ngrams(fasta)
