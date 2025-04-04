#!/usr/bin/env python3


seq="AGAATTGCGGAACTCTCGTAGTATTCTCATCTCACTCACCCGTCGCTACAGGCTTACGTGATGAGTCTAGTATACCAATCGAGTCTTCATTTTGATTTTAGAAATAGTTCCTACCCTTTCGCTTTTTCCAGTAGTTACCAGTTCAGGGCCGATGAGCCAGAGTCCTCAGACGCGATATTCACCCTACGGCTATCCAAACATTTGCGCTGTGCATTTGACGATCCGTTCAAGTAAAGCGAGCGGATGCGTGCGCAAGTACCCTTTTGCCCAAGACGATCAAGTGTCGATTGTTGGTGTCTCACGCCAGCAAGGTTGGGCGGACATGACCACGGACCAGTAAAAGTTACAATTCTATGGACATAGTATCTAAAGTGGCTTTTCAATGGAGATTAAGGGCCAAATGAGCAGCAGACATCTCGAACGCTTTCGAGGTTGTAACCCTTCCCGCTGCCGTTCCCTTATCTGTTATGAACCACATCTATTGCGCATTAAACGTTCCAACTAACAGTACGCGAAAGGATACCTCAGCGAGCTCCCCCGAATCGTGTTGGGGTGGTTTCTACACGGACTACATTTTAGCGAGAAGTTGCAATCCGTGCTATGCAGAAGCAATGCATCTCGGGGCGCGAGAGTGCACAACTAGCAGCAGTACACATTAGGCCGGCGGTAGCGCGATCAACGGTGAAATAGGGGCGTCCTCATTTTTACTTCTCCACACGAATAAATTCTTCGGGGAATAACGAAAAAGGATGCTCACCAGCATACGTGCCACTGTAGAGCAGCTAGGATGCCCCCGTCCGGGCACCGACCTAGGGTGGTTGCCCCAAAGTTACATGACAAATTTTATCCCTGCGGTTACAACCGGAGAGTGCAACAAAGCACGTTGAGTGGCCTCCCCTGAG"

def main():
    count = {}
    for i in seq:
        count[i] = count.get(i, 0) + 1
    print(count)

if __name__ == "__main__":
    main()
