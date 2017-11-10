#! /usr/bin/env python3

class Sequence:
    # Detect the type of sequence this is
    def detect_type(seq):
        # Check to see if this seq is only ACGT
        count = seq.count("A") + seq.count("C") + seq.count("G") + seq.count("T")

        return "N" if count == len(seq) else "P"

    def get_aa():
    # Obtained from https://stackoverflow.com/questions/19521905/translation-dna-to-protein
        return {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
        }

    # Perform the reverse complement on a nucleotide sequence
    def rev_comp(seq):
        ret_seq = ""
        rev_lookup = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

        for nuc in seq[::-1]:
            ret_seq += rev_lookup.get(nuc, '')

        return ret_seq

    # Translate the sequence on the specified frame
    def translate_frame(seq, frame):
        ret_protein = ""
        aa_lookup = Sequence.get_aa()

        # Frames 0-2 fwd, 3-5 rev
        strand = seq if frame < 3 else Sequence.rev_comp(seq)

        # Start on the requested frame, and skip a codons length each time
        for pos in range(frame % 3, len(strand), 3):
            codon = strand[pos:pos+3]

            if codon not in aa_lookup:
                ret_protein += "X"
            else:
                ret_protein += aa_lookup[codon]

        return ret_protein

    def __init__(self, seq = '', molecule = 'N'):
        # Set sequence (no override necessary)
        self.seq = seq
        self.molecule = molecule

        # Perform a type check if nucleotide
        if self.molecule == 'N':
            self.molecule = Sequence.detect_type(self.seq)

    # String method returns format "<TYPE>:<SEQ>"
    def __str__(self):
       return "{0}:{1}".format(self.molecule, self.seq)

    # Len should return length of sequence
    def __len__(self):
        return len(self.seq)

    # Equality checks for same sequence (or translated sequence for N==P)
    def __eq__(self, operand):
        # If not a sequence, return false
        if not isinstance(operand, Sequence):
            return False

        if self.molecule == operand.molecule:
            # If the same type, simply compare strings
            return self.seq == operand.seq
        else:
            # If different types, compare translations
            for trans1 in self.translate():
                for trans2 in operand.translate():
                    if trans1 == trans2:
                        return True

        return False

    # Addition should create a new sequence that's the concatenation of two
    # sequences of same type, or an empty nucleotide sequence if not
    def __add__(self, operand):
        # If not a sequence, return new empty sequence
        if not isinstance(operand, Sequence):
            return Sequence()

        # If not same type, return new empty sequence
        if self.molecule != operand.molecule:
            return Sequence()

        # Return concatenation of same type
        return Sequence(self.seq + operand.seq, self.molecule)

    # If nucleotides, translate on all 6 frames.  If prot, return
    def translate(self):
        ret_list = []

        if self.molecule == 'P':
            # For proteins, simply return sequence
            ret_list.append(self.seq)
        else:
            # For nucleotides, translate all 6 frames and add to list
            for frame in range(6):
                ret_list.append(Sequence.translate_frame(self.seq, frame))

        return ret_list

    # Calculate GC content
    def gc_content(self):
        # For proteins, return 0
        if self.molecule == 'P':
            return 0

        return (self.seq.count("C") + self.seq.count("G")) / len(self)

# Prompt the user for a list of nucleotide sequences
def get_seqs():
    seq_list = []

    while True:
        seq = input("Enter nucleotide seq: ")
        if seq == "":
            break

        seq_list.append(Sequence(seq))

    return seq_list

# Count kmers across all Sequence instances in the list
def get_kmers(seqs, size):
    kmers = {}

    for seq in seqs:
        # Skip protein sequences
        if seq.molecule == 'P': continue

        # Iterate through the starting position of each kmer and record
        for pos in range(len(seq) - size + 1):
            kmer = seq.seq[pos:pos+size]
            kmers.setdefault(kmer, 0)
            kmers[kmer] += 1

    return kmers

# Get seqs, and print 3mers, 4mers and 5mers
seqs = get_seqs()

print(get_kmers(seqs, 3))
print(get_kmers(seqs, 4))
print(get_kmers(seqs, 5))
