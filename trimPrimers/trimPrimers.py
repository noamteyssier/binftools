#!/usr/bin/env python3

import sys, os, argparse, gzip, re
from datetime import datetime
import multiprocessing as mp
import subprocess as sp
import triegex


class TrimPrimers:
    def __init__(self, parser):
        parser.parse_args(namespace=self) # add args to class

        self.gen_r1 = None
        self.gen_r2 = None
        self.gen_p1 = None
        self.gen_p2 = None
        self.regex_p1 = None
        self.regex_p2 = None
        self.trim_1 = None
        self.trim_2 = None

        self.seq_memory = dict()
        self.trim_memory = dict()
    def seq_gen(self, fname, type='fastq'):
        """sequence generator for fastq and fasta"""
        f_gen = self.open_file(fname)
        if type=='fastq':
            while True:
                try:
                    header, seq, _, quality = [next(f_gen).strip() for _ in range(4)]
                    yield [header, seq, quality]
                except StopIteration:
                    break

        elif type=='fasta':
            while True:
                try:
                    header, seq = [next(f_gen).strip() for _ in range(2)]
                    yield [header, seq]
                except StopIteration:
                    break
    def open_file(self, fname, how='read'):
        """return open file object for for gzip or non"""
        open_method = 'r' if how == 'read' else 'w'
        if fname[-3:] == '.gz':
            return gzip.open(fname, open_method + 't')
        else:
            return open(fname, open_method)
    def construct_primer_regex(self):
        """create regex expression for primer sets"""
        self.gen_p1, self.gen_p2 = [self.seq_gen(i, type='fasta') for i in [self.primer1, self.primer2]]
        self.regex_p1 = triegex.Triegex(list('^' + seq for (header, seq) in self.gen_p1)).to_regex()
        self.regex_p2 = triegex.Triegex(list('^' + seq for (header, seq) in self.gen_p2)).to_regex()
    def apply_regex(self, regex, sequence):
        """apply regex on sequence if not done before"""
        # only apply regex if not done before
        if sequence not in self.seq_memory:
            match = re.match(regex, sequence)
            self.seq_memory[sequence] = match

        return self.seq_memory[sequence]
    def apply_trim(self, record, match):
        """apply trim on fastq sequence and quality scores if not done before"""
        if record[1] not in self.trim_memory:
            start, end = match.span()
            self.trim_memory[record[1]] = [record[0], record[1][end:], record[2][end:]]
        return self.trim_memory[record[1]]
    def write_fastq(self, f, record):
        """write record in fastq format"""
        f.write("{0}\n{1}\n+\n{2}\n".format(record[0], record[1], record[2]))
    def Trim(self):
        """
        main method of class
        - begins multiprocessing reader and writer
        - initialize sequence generators
        - initialize regex expressions
        - start main loop
            - generator + regex
            - write
        """

        self.gen_r1, self.gen_r2 = [self.seq_gen(i) for i in [self.read1, self.read2]]
        self.construct_primer_regex()

        # standard
        trim_f1 = self.open_file(self.trim1, how='write')#.replace('.gz', ''), how='write')
        trim_f2 = self.open_file(self.trim2, how='write')#.replace('.gz', ''), how='write')

        while True:
            try:
                record1, record2 = [next(g) for g in [self.gen_r1, self.gen_r2]]
                match1 = self.apply_regex(self.regex_p1, record1[1])
                match2 = self.apply_regex(self.regex_p2, record2[1])

                if match1 and match2:
                    t1, t2 = [self.apply_trim(r, m) for (r,m) in [(record1, match1), (record2, match2)]]
                    self.write_fastq(trim_f1, t1)
                    self.write_fastq(trim_f2, t2)

            except StopIteration:
                break

def get_args():
    p = argparse.ArgumentParser()
    p.add_argument('-r1', '--read1', required=True)
    p.add_argument('-r2', '--read2', required=True)
    p.add_argument('-p1', '--primer1', required=True)
    p.add_argument('-p2', '--primer2', required=True)
    p.add_argument('-t1', '--trim1', required=True)
    p.add_argument('-t2', '--trim2', required=True)
    return p
def main():
    parser = get_args()
    t = TrimPrimers(parser)
    t.Trim()




    pass



if __name__ == '__main__':
    main()
