#!/usr/bin/env python3



def transcribe(seq):
    tmp=""
    for i in seq:
        if i == "T":
            tmp += "U"
        else:
            tmp += i
    return tmp
def main(seq):
    rna = transcribe(seq)
    return(rna)

if __name__ == "__main__":
    main()
