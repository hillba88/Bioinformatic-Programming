#! /usr/bin/env python3

STOP_CODONS = ['TAG', 'TAA', 'TGA']

# Reverse the string, then complement each nucleotide
def rev_comp(seq):
    ret_seq = ""

    for nuc in seq[::-1]:
        if nuc == 'A':
            ret_seq += 'T'
        if nuc == 'C':
            ret_seq += 'G'
        if nuc == 'G':
            ret_seq += 'C'
        if nuc == 'T':
            ret_seq += 'A'

    return ret_seq

# Check for stop codons
def check_stop(seq):
    # Check for each codon in forward sequence
    for codon in STOP_CODONS:
        if codon in seq:
            return True

    # Check for each codon in reverse sequence
    reverse = rev_comp(seq)

    for codon in STOP_CODONS:
        if codon in reverse:
            return True

    return False

# Calculate GC content
def gc_content(seq):
    gc_count = seq.count("G") + seq.count("C")
    return gc_count / len(seq)

# Master function, call one of the other functions depending on command
def master(seq, command):
    result = None

    # Reverse complement
    if command == 'revcomp':
        result = rev_comp(seq)

    # Stop codon detection
    if command == 'stop':
        result = check_stop(seq)

    # GC Content
    if command == 'gc':
        result = gc_content(seq)

    return result

seq = input("Enter a sequence: ")

# Try each command
print("Reverse Complement: ", master(seq, 'revcomp'))
print("Stop Codon Detect: ", master(seq, 'stop'))
print("GC Content: ", master(seq, 'gc'))
print("Invalid: ", master(seq, 'bad!'))
