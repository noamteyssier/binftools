#!/usr/bin/python



"""
A processing script to transform snp positions to regional snpStrings

input (tab delim)   : [chr pos base window strain]
output (tab delim)  : [strain window snpString]

"""

import sys, argparse

def parseTab(t):
    """yield parsed tab delim file"""
    with open(t) as f:
        next(f) # skip header
        while True:
            try:
                yield next(f).strip('\n').split('\t')
            except StopIteration:
                break
def catSNP(t):
    d = dict() # strain : window : [base i..j]
    for line in parseTab(t):
        chr, pos, base, window, strain = line
        if strain not in d:
            d[strain] = dict()
        if window not in d[strain]:
            d[strain][window] = list()
        d[strain][window].append(base)

    return d

def printout(d):
    for strain in d:
        for window in d[strain]:
            items = [strain, window, ''.join(d[strain][window])]
            print '\t'.join(items)




def main():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--tab', help = 'tab delim [ chr pos base window strain ]', required = True)
    args = p.parse_args()
    d = catSNP(args.tab)
    printout(d)


if __name__ == '__main__':
    main()
