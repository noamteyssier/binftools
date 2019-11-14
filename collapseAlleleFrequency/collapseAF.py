#!/usr/bin/env python3

import numpy as np
from tqdm import tqdm
import gzip


class CollapseAF(object):
    """

    Collapse stacked allele frequency tables based with non unique
    genomics positions

    Parameters
    ----------
    **kwargs : type
        keyword specific updates

    Attributes
    ----------
    skip_header : type
        skips first line of file
    sites : type
        dictionary to hold positions and arrays associated
    """

    def __init__(self, **kwargs):
        self.skip_header = True
        self.sites = {}

        self.__dict__.update(**kwargs)

    def reader(self, fn):
        """
        Generator of lines over a file

        Parameters
        ----------
        fn : type
            filename to iterate

        Yields
        -------
        list
            row split by tabs

        """

        def open_file(fn):

            if 'gz' in fn:
                return gzip.open(fn, 'rt')

            else:
                return open(fn, 'r')

        count = 0
        for line in open_file(fn):

            if self.skip_header and count == 0:
                count += 1
                continue

            if 'depth' in line:
                continue

            yield line.strip('\n').split('\t')

    def collapse(self, fn):
        """
        Collapses allele frequencies

        Parameters
        ----------
        fn : str
            filename to collapse

        Returns
        -------
        method
            results method

        """

        for line in tqdm(self.reader(fn)):

            site = line[0]
            depth = int(line[1])

            if depth > 0:
                frequencies = np.array(line[2:]).astype(float)
            else:
                frequencies = np.zeros(4)

            if site not in self.sites:
                self.sites[site] = np.zeros(4)

            self.sites[site] += (depth * frequencies)

        self.results()

    def results(self):
        """
        print results to stdout
        """

        for site in self.sites:
            depth = self.sites[site].sum()
            freqs = self.sites[site] / depth

            row = [site, depth] + freqs.tolist()

            print(
                '\t'.join([str(r) for r in row])
            )


def main():
    fn = 'fullset.allelefreq.tab.gz'
    c = CollapseAF()
    c.collapse(fn)


if __name__ == '__main__':
    main()
