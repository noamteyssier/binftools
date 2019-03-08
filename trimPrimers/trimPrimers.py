#!/usr/bin/env python3

import sys, os, argparse, gzip, re
import multiprocessing as mp
import subprocess as sp


class TrimPrimers:
    def __init__(self, parser):
        parser.parse_args(namespace=self) # add args to class

        self.gen_r1 = None
        self.gen_r2 = None
        self.gen_p1 = None
        self.gen_p2 = None
        self.regex_p1 = None
        self.regex_p2 = None
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
        self.regex_p1, self.regex_p2 = ['|'.join(['^' + seq for (header, seq) in g])
            for g in [self.gen_p1, self.gen_p2]]
    def apply_trim(self, record, match):
        """apply trim on fastq sequence and quality scores"""
        start, end = match.span()
        return [record[0], record[1][end:], record[2][end:]]
    def write_fastq(self, fname, q):
        f = self.open_file(fname.replace('.gz',''), how='write')
        while True:
            record = q.get()
            if record != 'exit':
                f.write(
                    "{0}\n{1}\n+\n{2}\n".format(record[0], record[1], record[2])
                )
            else:
                break
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
        # mp.set_start_method('spawn')
        q1 = mp.Queue()
        q2 = mp.Queue()

        trim1_writer = mp.Process(target=self.write_fastq, args=(self.trim1, q1))
        trim2_writer = mp.Process(target=self.write_fastq, args=(self.trim2, q2))
        trim1_writer.start()
        trim2_writer.start()

        self.gen_r1, self.gen_r2 = [self.seq_gen(i) for i in [self.read1, self.read2]]
        self.construct_primer_regex()
        while True:
            try:
                record1, record2 = [next(g) for g in [self.gen_r1, self.gen_r2]]
                match1, match2 = [re.match(primers, seq) for (primers, seq) in [(self.regex_p1, record1[1]), (self.regex_p2, record2[1])]]

                if match1 and match2:
                    t1, t2 = [self.apply_trim(r, m) for (r,m) in [(record1, match1), (record2, match2)]]
                    q1.put(t1)
                    q2.put(t2)

            except StopIteration:
                q1.put('exit')
                q2.put('exit')
                trim1_writer.join()
                trim2_writer.join()
                if self.trim1[-3:] == '.gz':
                    sp.Popen(
                        "gzip -f {0} {1}".format(self.trim1.replace('.gz', ''), self.trim2.replace('.gz', '')),
                        shell=True
                        )
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
