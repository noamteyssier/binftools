#!/usr/bin/env python3

import numpy as np
import pandas as pd
import pysam as ps
import argparse, sys, os, shutil, gzip

class RegionExtractor:
    def __init__(self, bam_fn, bed_fn, output_dir, overwrite=False):
        self.input_bam = bam_fn
        self.input_bed = bed_fn

        self.bam = None
        self.bed = None
        self.output_dir = output_dir
        self.overwrite = overwrite

        self.intervalStats = {
            'interval_name' : [],
            'n_partial' : [],
            'n_total' : []
        }

        self._load_bam()
        self._load_bed()
        self._prepare_output_dir()
    def _eprint(self, *args, **kwargs):
        """stderr from https://stackoverflow.com/a/14981125/8767800"""
        if self.verbose:
            print(*args, file=sys.stderr, **kwargs)
    def _load_bam(self):
        """load bam into class"""
        self.bam = ps.AlignmentFile(self.input_bam, "rb")
    def _load_bed(self):
        """load region bed into class"""
        self.bed = pd.read_csv(self.input_bed, sep="\t", header=None)
        self.original_bed_size = self.bed.shape[1]
        self.bed['region'] = self.bed[0] + ":" + self.bed[1].astype('str') + '-' + self.bed[2].astype('str')
    def make_directory(self, name):
        """
        Check if directory exists
        Check if overwrite flag
        Create directory if flag given or nonexistant
        """

        if os.path.exists(name):
            if not self.overwrite:
                sys.exit(
                    'ERROR  : directory already exists : {}\nrerun with --overwrite'.\
                    format(name)
                    )
            else:
                shutil.rmtree(name)
        os.makedirs(name)
    def _prepare_output_dir(self):
        """
        create default output directory name if necessary
        create directory name
        """
        if not self.output_dir:
            self.output_dir = self.input_bam.replace('.bam', '_out/')
        if self.output_dir[-1] != '/':
            self.output_dir = self.output_dir + '/'
        self.make_directory(self.output_dir)
    def _call_pileup(self, interval):
        """call pileup given interval object"""
        return self.bam.pileup(
            contig=interval[0],
            start=interval[1],
            stop=interval[2],
            truncate=True
            )
    def getIntervalOFN(self, interval):
        """
        prepare interval name if not given in bed and interval filename
        """
        try:
            interval_name = str(interval[3])
        except KeyError:
            interval_name = '.'.join([str(i) for i in interval[:3]])
        interval_ofn = self.output_dir + str(interval_name) + '.fastq.gz'
        return interval_name, interval_ofn
    def GetIntervalInfo(self, interval):
        """
        organize name, sequence, and quality as tuple for each interval base position
        """
        return [(
            np.array(column.get_query_names()),
            np.array(column.get_query_sequences()),
            np.array(column.get_query_qualities())
            ) for column in self._call_pileup(interval)]
    def UpdateStatistics(self, interval_name, total_headers, passing_headers):
        """
        update dictionary entry for interval name
        """
        self.intervalStats['interval_name'].append(interval_name)
        self.intervalStats['n_partial'].append(total_headers.size)
        self.intervalStats['n_total'].append(passing_headers.size)
    def SelectFullSpanning(self, tuples):
        """
        Select header names that are present across all interval positions
        """
        if len(tuples) == 0:
            return [np.array([])] * 2

        headers = np.hstack([t[0] for t in tuples])
        unique_headers, base_coverage = np.unique(headers, return_counts = True)
        passing_headers = unique_headers[np.where(base_coverage == len(tuples))[0]]
        return unique_headers, passing_headers
    def BuildSequences(self, tuples, passing_headers):
        """
        - Iterate across columns
        - Select indices with fully spanning reads
        - return dictionary of header : sequence, quality
        """
        read_dictionary = {h : {'seq' : [], 'qual' : []} for h in passing_headers}
        for t in tuples:
            idx = np.isin(t[0], passing_headers)
            for i, header in enumerate(t[0][idx]):
                read_dictionary[header]['seq'].append(t[1][i].upper())
                read_dictionary[header]['qual'].append(t[2][i])
        return read_dictionary
    def ParseReadDictionary(self, read_dictionary):
        """
        parse read dictionary and return header, sequence, and quality as list
        """
        for i in read_dictionary:
            header = i
            seq = ''.join(read_dictionary[i]['seq'])
            qual = ps.qualities_to_qualitystring(read_dictionary[i]['qual'])
            yield [header, seq, qual]
    def WriteSequences(self, read_dictionary, interval_ofn):
        """open interval file and write fastq format"""
        interval_file = gzip.open(interval_ofn, 'wt')
        for header, seq, qual in self.ParseReadDictionary(read_dictionary):
            interval_file.write(
                '@{}\n{}\n+\n{}\n'.\
                format(header, seq, qual)
                )
    def GetInterval(self, interval):
        """
        Extract interval from bam file
        """
        interval_name, interval_ofn = self.getIntervalOFN(interval)
        tuples = self.GetIntervalInfo(interval)
        total_headers, passing_headers = self.SelectFullSpanning(tuples)
        self.UpdateStatistics(interval_name, total_headers, passing_headers)

        # only write file if reads found
        if passing_headers.size > 0:
            read_dictionary = self.BuildSequences(tuples, passing_headers)
            self.WriteSequences(read_dictionary, interval_ofn)
    def PrintStats(self):
        interval_stats = pd.DataFrame.from_dict(self.intervalStats)
        interval_stats.to_csv(sys.stdout, sep="\t", index=False)
        interval_stats.to_csv(self.output_dir + 'intervalStats.tab', sep="\t", index=False)
    def Extract(self):
        """iterate over intervals and run extraction"""
        self.bed.apply(
            lambda x : self.GetInterval(x), axis=1
        )
        self.PrintStats()



def get_args():
    p = argparse.ArgumentParser()
    p.add_argument('-i',
        '--input_bam',
        help='input bam to extract regions',
        required=True
        )
    p.add_argument('-b',
        '--input_bed',
        help='input bed to extract regions',
        required=True
        )
    p.add_argument('-o',
        '--output_dir',
        help='output directory name to write to (default=$BAM_out)',
        required=False,
        default=None
        )
    p.add_argument('--overwrite',
        help='overwrite output directory',
        required=False,
        action='store_true')
    args = p.parse_args()
    return args
def main(args):
    re = RegionExtractor(args.input_bam, args.input_bed, args.output_dir, args.overwrite)
    re.Extract()

if __name__ == '__main__':
    args = get_args()
    main(args)
