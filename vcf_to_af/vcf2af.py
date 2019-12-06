#!/usr/bin/env python3

import numpy as np
import pandas as pd
import pysam
import argparse
import sys
from tqdm import tqdm

class vcf2af(object):
    """
    Class to convert variant call files to allele frequency tables

    requires :
        - variant call file
        - interval bed file
    """

    def __init__(self, vcf_fn, bed_fn):
        self.vcf_fn = vcf_fn
        self.bed_fn = bed_fn

        self.init_vcf()

    def init_vcf(self):
        """
        load vcf into class, isolates sample names from header
        """

        self.vcf = pysam.VariantFile(self.vcf_fn)
        self.samples = [s for s in self.vcf.header.samples]

    def bed_iter(self):
        """
        generator of intervals from bed file
        """

        with open(self.bed_fn) as f:
            while True:
                try:
                    line = next(f).strip('\n')
                    yield line.split('\t')
                except StopIteration:
                    break

    def get_allele_frequencies(self, pos):
        """
        create allele frequencies from allelic depths for a given position
        """

        allele_depths = [
            np.array(s['AD']) for s in pos.samples.values()
            ]

        allele_freqs = [
            ad/ad.sum() for ad in allele_depths
            ]

        return allele_freqs

    def process_position(self, pos, position_name):
        """
        create dataframe of allele frequencies for a position for each sample
        """

        allele_freqs = self.get_allele_frequencies(pos)

        num_alleles = np.max(
            [af.size for af in allele_freqs]
            )

        pos_frame = pd.DataFrame(
            np.vstack(allele_freqs)
            )

        pos_frame['sample'] = self.samples
        pos_frame['locus'] = position_name

        pos_frame = pos_frame.melt(
            id_vars=['locus', 'sample'],
            var_name = 'allele',
            value_name = 'fraction'
            )

        pos_frame = pos_frame[pos_frame.fraction > 0]

        return pos_frame

    def process_missing(self, position_name):
        """
        create a dataframe for missing intervals for each sample
        with a default reference allele
        """

        frame = pd.DataFrame([
            [position_name] * len(self.samples),
            self.samples,
            [0] * len(self.samples),
            [1.0] * len(self.samples)
            ]).T
        frame.columns = ['locus', 'sample', 'allele', 'fraction']
        return frame

    def fit(self):
        """
        creates a dataframe of
        locus, sample, allele, fraction
        """

        allele_frame = []

        for b in self.bed_iter():
            chrom, start, end = b[:3]
            position_name = '.'.join(b)

            region_call = self.vcf.fetch(chrom, int(start), int(end))

            try:
                pos = next(region_call)
                allele_frame.append(
                    self.process_position(pos, position_name)
                    )

            except StopIteration:
                # missing region

                allele_frame.append(
                    self.process_missing(position_name)
                )

        allele_frame = pd.concat(allele_frame)

        allele_frame.to_csv(sys.stdout, sep="\t", index=False)

def main(args):

    v = vcf2af(
        vcf_fn = args.input_vcf,
        bed_fn = args.input_bed
        )
    v.fit()

    pass

def get_args():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input_vcf',
        required=True,
        help='input vcf to calculate AF'
        )
    p.add_argument('-b', '--input_bed',
        required=True,
        help='input bed to intersect'
        )

    args = p.parse_args()
    return args

if __name__ == '__main__':
    args = get_args()
    main(args)
