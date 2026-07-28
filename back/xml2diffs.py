#!/usr/bin/env python

import sys
from Bio.Blast import NCBIXML

def update_diffs(diffs, conserv, hsp):
    pos = -1
    while True:
        pos = hsp.match.find(conserv, pos + 1)
        if pos < 0:
            break
        diffs.append(pos)


def construct_diffs(parsed):
    diffs = []
    for record in parsed:
        if record.alignments:
            count_align = 0
            for align in record.alignments:
                count_align += 1
                if count_align > 1:
                    print('More than one hit for ' + record.query.split(' ')[1], file=sys.stderr)
                    break
                count_hsp = 0
                for hsp in align.hsps:
                    with open('output', 'w') as f:
                        print(hsp.match, file=f)
                        print(hsp.sbjct, file=f)
                    count_hsp += 1
                    if count_hsp > 1:
                        print('More than one hsp for ' + record.query, file=sys.stderr)
                    update_diffs(diffs, ' ', hsp)
                    update_diffs(diffs, '+', hsp)
    return diffs

def diffs(xml):
    parsed = NCBIXML.parse(open(xml))

    return construct_diffs(parsed)

                            
def main():
    XML = sys.argv[1]

    print(diffs(XML))


if __name__ == "__main__":
    main()
