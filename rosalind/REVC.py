#!/usr/bin/env python3

#provide complement dictionary
comp_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}


def rev(seq):
    return seq[::-1]


def complement(seq):
    ## to find seq complement
    comp = ''
    for i in seq:
        comp += comp_dict[i]
    return comp

def main(seq): 
    seq = rev(seq)
    seq = complement(seq)
    return seq

if __name__ == "__main__":
    main()
