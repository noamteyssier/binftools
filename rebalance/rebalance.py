#!/usr/bin/env python3

import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt
import argparse, sys
sns.set(rc={'figure.figsize':(22, 22), 'lines.linewidth': 5})
pd.options.mode.chained_assignment = None  # default='warn'

class RebalancePool(object):
    """
    Filter a dataset based on minimum reads and fractions
    Rebalance the filtered pool based on mapping fraction
    """
    def __init__(self, flagstat_frame, minimum_reads=6e4, minimum_fraction=0.2, scale=10, plot = False):
        self.frame = pd.read_csv(flagstat_frame, sep="\t", header=None)
        self.frame['fraction'] = self.frame.iloc[:,1] / self.frame.iloc[:,0]

        self.minimum_reads = minimum_reads
        self.minimum_fraction = minimum_fraction
        self.scale = scale
        self.plot = plot

        if self.plot:
            self.plot_distributions()

    def distribution_totalReads(self):
        """
        Plot distribution of total reads across all samples
        """
        sns.distplot(
            np.log10(self.frame.iloc[:,0])
            )
        plt.axvline(
            np.log10(self.minimum_reads),
            ls = ':', color = 'black'
        )
        plt.show()
    def distribution_readFraction(self):
        """
        Plot distribution of read fractions across all samples
        """
        sns.distplot(
            self.frame.fraction
            )
        plt.axvline(
            self.minimum_fraction,
            ls = ':', color = 'black'
        )
        plt.show()
    def plot_distributions(self):
        """
        Plot distributions and overlay filters
        """

        self.distribution_totalReads()

        self.distribution_readFraction()
    def filter(self):
        """
        Apply read count and fraction filters
        """
        self.filtered = self.frame[
            (self.frame.iloc[:,0] >= self.minimum_reads) & \
            (self.frame.fraction >= self.minimum_fraction)
            ]
        return self.filtered
    def fit(self):
        """
        Apply rebalancing across pool after filtering
        """
        self.filter()
        self.filtered['rebalance'] = (
            1 / (self.filtered.fraction / self.filtered.fraction.sum()) / self.scale
        )
        return self.filtered


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("-i", '--input_file',
        help='Input flagstat to process (total, mapped, paired, sampleName)',
        required=True
        )
    p.add_argument("-r", '--minimum_reads',
        help='minimum read threshold',
        required=False,
        default = 6e4,
        type=float
        )
    p.add_argument("-f", '--minimum_fraction',
        help='minimum mapping fraction threshold',
        required=False,
        default = 0.2,
        type=float
        )
    p.add_argument("-s", '--scale',
        help='scalar to put rebalance with',
        required=False,
        default = 10,
        type=float
        )
    p.add_argument('-p', '--plot',
        help='plot thresholds and distributions',
        required=False,
        action='store_true'
        )
    args = p.parse_args()
    return args
def main(args):
    rp = RebalancePool(
        args.input_file,
        minimum_reads = args.minimum_reads,
        minimum_fraction = args.minimum_fraction,
        scale = args.scale,
        plot = args.plot
    )


    rebalance = rp.fit()
    rebalance.to_csv(sys.stdout, sep="\t", index=False)

if __name__ == '__main__':
    args = get_args()
    main(args)
