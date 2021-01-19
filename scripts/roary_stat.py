#/usr/bin/python3.7

import sys
import os
import csv

constant_cols = 14

def main():
    filename=sys.argv[1]

    with open(filename) as gpa:
        cr = csv.reader(gpa)
        headers = next(cr, None)
        if headers:
            nois = len(headers) - constant_cols
            stats = [[0 for x in range(nois)] for x in range(nois)]
            names = headers[constant_cols:]
        else:
            print('Empty file', file=sys.stderr)
            sys.exit()

        for row in cr:
            n = int(row[3])
            for i in range(nois):
                if row[constant_cols + i]:
                    stats[n-1][i] += 1

    header = ['Count'] + names
    print('\t'.join(header))
    for i, row in enumerate(stats):
        rowstr = [str(i+1)] + [str(x) for x in row]
        print('\t'.join(rowstr))


if __name__ == '__main__' :
    main()
