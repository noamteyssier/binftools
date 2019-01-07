#!/usr/bin/env python3

import sys
import argparse

def get_args():
    """method to handle arguments"""
    p = argparse.ArgumentParser()
    p.add_argument('-o', '--output_file', help = 'output file to write to (default = STDOUT)')
    p.add_argument('-n', '--name', help = 'name to to append as column')
    args = p.parse_args()
    return args
def print_out(values, output_fn=None):
    """write to file or stdout"""
    line = '\t'.join(values)
    if output_fn:
        open(output_fn, 'w+').write(line + '\n')
    else:
        print(line)
def main():
    args = get_args()
    relevant_info = []
    for line in sys.stdin:
        if 'total' in line:
            relevant_info.append(line.split(' ')[0])
            continue
        elif 'mapped' in line and 'N/A' in line:
            relevant_info.append(line.split(' ')[0])
            continue
        elif 'properly paired' in line:
            relevant_info.append(line.split(' ')[0])
            continue

    if args.name:
        relevant_info.append(args.name)

    print_out(relevant_info, args.output_file)



if __name__ == '__main__':
    main()
