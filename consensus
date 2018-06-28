#!/usr/bin/python
import sys, argparse


class Region:
    """class to handle regions and store bases"""
    def __init__(self, entry):
        if len(entry) == 5:
            self.chrom, self.left, self.right, self.name, self.direction = entry
        else:
            try:
                self.chrom, self.left, self.right, self.name = entry
            except ValueError:
                sys.exit("ERROR : bed format must have at least 4 values (chr, left, right, label) and up to 5")
        self.entry = entry
        self.bases = []
    def assignPos(self, line, threshold):
        """assign a base to class if it is within chromosomal region"""
        ch,pos,ref,depth,pile,cig = line
        if self.chrom == ch:
            if (pos >= self.left) and (pos <= self.right):
                self.bases.append(parse_pile(pile, depth, threshold))
    def printConsensus(self, label):
        """print formatted output"""
        items = self.entry + [''.join(self.bases)]
        #items = [self.chrom, self.left, self.right, self.name, self.direction, ''.join(self.bases)]
        if label:
            items.append(label)
        print '\t'.join(items)

def max_hist(string, lex=False):
    """take histogram of each character in string and return highest"""
    hist = dict()
    for s in string:
        if lex: # add lexicon to check against
            if s not in lex:
                continue
        if s not in hist:
            hist[s] = 0
        hist[s] += 1
    return max(hist, key = hist.get)
def parse_stdin():
    """generate stripped tab delim lines"""
    for line in sys.stdin:
        yield line.strip('\n').split('\t')
def parse_pile(pile, depth, threshold):
    """turn pile into uppercase and return maximum"""
    pile = pile.upper()
    lex = ['A', 'C', 'G', 'T']
    if int(depth) >= int(threshold):
        return max_hist(pile, lex)
    else:
        return 'x'
def parse_bed(bed):
    with open(bed) as f:
        while True:
            try:
                entry = next(f).strip('\n').split('\t')
                yield Region(entry)
            except StopIteration:
                break
def makeBED(bed):
    return [r for r in parse_bed(bed)]
def writeout(ch, pos, pile, label=False):
    """formatted output"""
    items = [ch,pos,pile]
    if label:
        items.append(label)
    print '\t'.join(items)
def main():
    """
    Pull the consensus base from a pileup line by line
    output tab delim w/o bedArrange:
        chromosome position base
    output tab delim w/ bedArrange:
        chromosome leftBound rightBound regionName direction consensus
    """
    p = argparse.ArgumentParser()
    p.add_argument('-l', '--label', help='add label column to output')
    p.add_argument('-b', '--bedArrange', help='BED File to arrange bases')
    p.add_argument('-t', '--threshold', default = 100, help='minumum depth to generate consensus')
    args = p.parse_args()

    if sys.stdin.isatty(): # exit if no stdin
        sys.exit("ERROR: pass pileup through stdin")

    if not args.bedArrange:
        for line in parse_stdin():
            ch,pos,ref,depth,pile,cig = line
            writeout(ch, pos, parse_pile(pile, depth, args.threshold), args.label)
    else:
        regions = makeBED(args.bedArrange)
        for line in parse_stdin():
            [r.assignPos(line, args.threshold) for r in regions]
        [r.printConsensus(args.label) for r in regions]


if __name__ == '__main__':
    main()
