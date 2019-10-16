#!/usr/bin/env python3

import numpy as np
from tqdm import tqdm
from nltk import ngrams
from sklearn.decomposition import PCA
from sklearn.cluster import Birch, DBSCAN
from scipy.stats import bernoulli, mode

from multiprocess import Pool

import seaborn as sns
import matplotlib.pyplot as plt

import gzip
import sys


class EmbedSeq(object):

    def __init__(self, **kwargs):
        self.expected_length = None

        self.original_sequences = {}

        self.sequence_table = {}

        self.sequence_id = {}

        self.matrix = np.array([])

        self.alphabet = {
            'A': 0,
            'C': 1,
            'T': 2,
            'G': 3
            }

        self.index_alphabet = {
            0: 'A',
            1: 'C',
            2: 'T',
            3: 'G'
            }

        self.drop_n = True

        self.drop_length = True

        self.add_noise = True

        self.error_rate = 0.001

        self.verbose = True

        self.num_iter_exp_length = 300

        self.labels = np.array([])

        self.final_seqs = {}

        self.num_iterations = 100

    def AddNoise(self, seq):

        if seq not in self.original_sequences:
            self.original_sequences[seq] = 0

        self.original_sequences[seq] += 1

        noise = np.random.choice(
            list(self.alphabet.keys()),
            replace=True, size=self.expected_length
            )

        mask = bernoulli.rvs(
            self.error_rate, size=self.expected_length
            ).astype(bool)

        seq = np.array([i for i in seq])

        seq[mask] = noise[mask]

        return ''.join(seq)

    def iter_seq(self, fn):

        count = 0

        if self.verbose:
            iter = tqdm(self.sequence_reader(fn), desc='building embedding')
        else:
            iter = self.sequence_reader(fn)

        for record in iter:

            seq = record[1]

            if self.drop_length and len(seq) != self.expected_length:
                continue

            if self.drop_n and 'N' in seq:
                continue

            if self.add_noise:
                seq = self.AddNoise(seq)

            if seq not in self.sequence_table:
                self.sequence_table[seq] = 1
                self.sequence_id[count] = seq
                count += 1
                yield seq

            else:
                self.sequence_table[seq] += 1

    def fit(self, fn):

        if not self.expected_length:
            self.predict_expected_length(fn)

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

    def predict_expected_length(self, fn):

        lengths = np.zeros(self.num_iter_exp_length)

        for idx, rec in enumerate(self.sequence_reader(fn)):
            if idx == self.num_iter_exp_length:
                break
            lengths[idx] = len(rec[1])

        self.expected_length = int(mode(lengths).mode[0])

    def label(self, matrix, **kwargs):

        return self.__label__(matrix, **kwargs)


class seq_to_onehot(EmbedSeq):
    """
    Seq to One Hot Encoding
    """

    def __init__(self, **kwargs):

        EmbedSeq.__init__(self, **kwargs)

        self.flatten = True
        self.__dict__.update(kwargs)

    def __fit__(self, fn):

        self.matrix = np.stack(
            [self.convert(s) for s in self.iter_seq(fn)]
            )

    def __label__(self, matrix, **kwargs):

        c = ClusterEmbed(**kwargs)
        c.fit(matrix)
        c.label(method=Birch, n_clusters=None, plot=False)

        self.labels = c.labels
        return self.labels

    def convert(self, seq):

        arr = np.zeros(
            (self.expected_length, len(self.alphabet))
            )

        for i, base in enumerate(seq):
            arr[i][self.alphabet[base]] = 1

        if self.flatten:
            arr = arr.ravel()
        return arr

    def consensus(self, labels):

        def build_consensus():
            consensus = np.zeros(
                (np.unique(labels).size, self.expected_length, len(self.alphabet))
                )

            for idx, cluster in enumerate(labels):

                if self.flatten:
                    mat = self.matrix[idx].\
                        reshape(self.expected_length, len(self.alphabet))

                else:
                    mat = self.matrix[idx]

                consensus[cluster] += mat

            seqs = []
            for cluster in consensus:
                cluster_vals = cluster.argmax(axis=1)
                cluster_seq = ''.join([
                    self.index_alphabet[i] for i in cluster_vals
                    ])
                seqs.append(cluster_seq)

            return np.array(seqs)

        def validate_consensus_sequences(seqs):
            valid_seqs = np.array([
                s in self.original_sequences for s in seqs
                ])
            return valid_seqs

        consensus_seqs = build_consensus()
        valid_seqs = validate_consensus_sequences(consensus_seqs)

        wat = np.where(valid_seqs)[0]
        clusters, counts = np.unique(labels, return_counts=True)

        valid_counts = counts[wat]

        seqs, counts = np.unique(consensus_seqs[valid_seqs], return_counts=True)

        clustered_counts = np.zeros(len(seqs))

        for i, s in enumerate(seqs):
            idx = np.where(consensus_seqs[valid_seqs] == s)[0]
            c = valid_counts[idx].sum()
            clustered_counts[i] = c

        for s in seqs:
            if s not in self.final_seqs:
                self.final_seqs[s] = []

        for s in self.final_seqs:
            if s in seqs:
                self.final_seqs[s].append(
                    clustered_counts[np.where(seqs == s)[0]]
                    )
            else:
                self.final_seqs[s].append(np.zeros(1))

        return clustered_counts / clustered_counts.sum()


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


class ClusterEmbed(object):

    def __init__(self, **kwargs):
        self.n_components = 3

        self.pcs = np.array([])
        self.labels = np.array([])
        self.clusters = np.array([])
        self.fractions = np.array([])

    def __str__(self):
        return (
                "clusters : \n{}\nfractions : \n{}\n".format(
                    ' '.join([str(i) for i in self.clusters]),
                    ' '.join([str(i) for i in self.fractions])
                    )
                )

    def run_pca(self, mat):
        pca = PCA(n_components=self.n_components)
        pcs = pca.fit_transform(mat)
        return pcs

    def plot_pca(self, labels=None):

        sns.scatterplot(
            x=self.pcs[:, 0],
            y=self.pcs[:, 1],
            hue=labels
            )
        plt.show()

    def fit(self, mat):
        self.pcs = self.run_pca(mat)

    def transform(self):
        return self.pcs

    def fit_transform(self, mat):
        self.fit(mat)
        return self.transform()

    def label(self, method=Birch, plot=True, **kwargs):
        m = method(**kwargs)

        self.labels = m.fit_predict(self.pcs)

        self.clusters, self.counts = np.unique(
            self.labels, return_counts=True
        )
        self.fractions = self.counts / self.counts.sum()

        if plot:
            self.plot_pca(labels=self.labels)


def to_onehot(fasta):
    oh = seq_to_onehot(error_rate=0.04, num_iterations=5)
    mat = oh.fit_transform(fasta)

    for i in tqdm(range(10), desc='bootstrapping'):
        labels = oh.label(mat)
        fractions = oh.consensus(labels)
        print(fractions)
        print('')


    global_sum = 0
    for s in oh.final_seqs:
        oh.final_seqs[s] = np.concatenate(oh.final_seqs[s]).sum()
        global_sum += oh.final_seqs[s]

    for s in oh.final_seqs:
        print(
            oh.final_seqs[s] / global_sum
        )

    return fractions


def to_bag(fasta):
    ng = seq_to_bag(ngram_size=7, error_rate=0.01)
    mat = ng.fit_transform(fasta)

    c = ClusterEmbed(n_components=2)
    c.fit(mat)
    c.label(method=DBSCAN)
    print(c)


def run_tc5():
    fasta = "/home/noam/projects/embedCluster/seq_data/paired/TC5/t{}.fasta.gz"
    # target_arr = np.arange(1,101)
    # fasta_list = [fasta.format(i) for i in target_arr]

    target_arr = [1,2,3,10,76]
    fasta_list = [fasta.format(i) for i in target_arr]

    p = Pool()
    tc5 = p.map(to_onehot, fasta_list)

    for t, i in enumerate(tc5):
        print(
            '\t'.join(['t{}'.format(target_arr[t])] + [str(j) for j in i])
            )


if __name__ == '__main__':
    # run_tc5()
    # fasta = "/home/noam/projects/embedCluster/seq_data/paired/TC5/t{}.fasta.gz"

    fasta = 'downsampled_ts.fa'
    to_onehot(fasta)
