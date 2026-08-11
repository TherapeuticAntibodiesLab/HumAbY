#!/usr/bin/env python

import sys
from Bio import SeqIO
from Bio.Blast import NCBIXML

def update_diffs(diffs, conserv, hsp):
    pos = -1
    while True:
        pos = hsp.match.find(conserv, pos + 1)
        if pos < 0:
            break
        diffs.append(pos)


def diffsGerm(xml):
    parsed = NCBIXML.parse(open(xml))

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
                    mtch = hsp.match
                    sbjc = hsp.sbjct
                    count_hsp += 1
                    if count_hsp > 1:
                        print('More than one hsp for ' + record.query, file=sys.stderr)
                    update_diffs(diffs, ' ', hsp)
                    update_diffs(diffs, '+', hsp)
    return diffs, mtch, sbjc

def backmutate(aa_filename, opt1_dir, xml_filename, chain):
    print('\n', 'Chain', chain)
    
    ## Backmutations by Joao
    with open(aa_filename, 'r') as f:
        murine = f.readline().rstrip()
        if chain == 'VL':
            murine = f.readline().rstrip()

    ## get FR regions
    ## and print header: FR1-->  FR2-->  etc
    header = []
    fr_filename = opt1_dir + "/" + 'FRfile' + chain
    j = 0
    nfr = 1
    with open(fr_filename, 'r') as fr_file:
        fr_indices = []
        for line in fr_file:
            fr = line.rstrip()
            i = murine.find(fr)
            header.extend([' '] * (i - j))
            j = i + len(fr)
            label = 'FR' + str(nfr)
            nfr += 1
            header.extend(list(label))
            header.extend(['-'] * (len(fr) - len(label)))
            fr_indices.extend(list(range(i, i+len(fr))))
            FR4start = i
    print('      ', ''.join(header))

    ## get differences between murine and its main germline
    diffs, mtch, mGerm = diffsGerm(xml_filename)
    print('mGerm ', mGerm)
    print('match ', mtch)
    mur_name = aa_filename.split('/')[-1].split('.')[0][0:6]
    print(mur_name.ljust(6), murine)

    ## open human germline file with sequences
    hGerm_filename = opt1_dir + "/" + 'hGerm' + chain
    hGerm = open(hGerm_filename, 'r')
    ## loop for each proposed seq
    proposed_filename = opt1_dir + "/" + 'humanized_'
    proposed_filename += chain.lower() + '.fasta'
    for record in SeqIO.parse(proposed_filename, "fasta"):
        ## print human germline corresponding to proposed
        ## they are in the hGerm file, in thre same order
        line = hGerm.readline()
        print('hGerm ', line, end='')

        ## Now process backmutations & print result
        proposed = record.seq
        backmutated = list(proposed)
        applied_rules = {
            'murine_germline_difference': set(),
            'cysteine_proline': set(),
            'position_71': set(),
            'fr4_motif': set(),
        }

        ## for each nonzero difference in FR region:
        ##    proposed <- murine
        for i in diffs:
            if i > 0 and i in fr_indices:
                if backmutated[i] != murine[i]:
                    applied_rules['murine_germline_difference'].add(i + 1)
                backmutated[i] = murine[i]
        ## for each disagreement between murine and proposed in FR region:
        ##   if proposed == cysteine or proline:
        ##     proposed <- murine
        ##   if kabat == 71:
        ##     proposed <- murine
        for i in fr_indices:
            if proposed[i] != murine[i]:
                if proposed[i] == 'C' or proposed[i] == 'P':
                    if backmutated[i] != murine[i]:
                        applied_rules['cysteine_proline'].add(i + 1)
                    backmutated[i] = murine[i]
            if i == 71:
                if backmutated[i] != murine[i]:
                    applied_rules['position_71'].add(i + 1)
                backmutated[i] = murine[i]

        ## if proposed FR4 does not start with WGXG:
        ##    proposed FR4[0:1] <- WG
        ##    proposed FR4[3] <- G
        fr4_values = {
            FR4start: 'W' if chain == 'VH' else 'F',
            FR4start + 1: 'G',
            FR4start + 3: 'G',
        }
        for position, residue in fr4_values.items():
            if backmutated[position] != residue:
                applied_rules['fr4_motif'].add(position + 1)
            backmutated[position] = residue
        print('backmt', "".join(backmutated))
        for rule_name, positions in applied_rules.items():
            ordered = sorted(positions)
            print(
                'rule', rule_name,
                'count=' + str(len(ordered)),
                'positions=' + (','.join(map(str, ordered)) if ordered else '-')
            )
    hGerm.close()

def main():
    aa_filename = sys.argv[1]
    opt1_dir = sys.argv[2]
    xmlH_filename = sys.argv[3]
    backmutate(aa_filename, opt1_dir, xmlH_filename, 'VH')
    xmlL_filename = sys.argv[4]
    backmutate(aa_filename, opt1_dir, xmlL_filename, 'VL')

if __name__ == "__main__":
    main()
