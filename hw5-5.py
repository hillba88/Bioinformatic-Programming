#! /usr/bin/env python3


class Sequence:
    valid_dna = "ATCG"

    def get_lookup():
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

    # creates a helper class method which detects whether a given "seq" is a
    # nucleotide or protein sequence
    def sequence_type(seq):
        valid_seq = "ATCG"
        for i in seq:
            if i not in valid_seq:
                return False
        return True

    # returns reverse compliment to nucleotide sequence
    def rev_comp(self):
        ret_seq = ""

        for nuc in self.seq[::-1]:
            if nuc == 'A':
                ret_seq += 'T'
            if nuc == 'C':
                ret_seq += 'G'
            if nuc == 'G':
                ret_seq += 'C'
            if nuc == 'T':
                ret_seq += 'A'

        return ret_seq

    # translates a given nucleotide sequence into all possible amino acids
    # in all six reading frames
    def translateSeq(self):
        ret_protein = []
        lookup = Sequence.get_lookup()

        # Start on the requested frame, and then skip a codons length each time
        for frame in range(6):
            strand = self.seq if frame < 3 else self.rev_comp()
            trans_seq = ''

            for pos in range(frame % 3, len(strand), 3):
                if frame < 3:
                    codon = strand[pos:pos+3]
                    if codon in lookup:
                        trans_seq += lookup[codon]
                else:
                    rev_codon = strand[pos:pos+3]
                    if rev_codon in lookup:
                        trans_seq += lookup[rev_codon]

            ret_protein.append(trans_seq)


        # if protein sequence, return sequence. If nucleotide sequence, return translation on 6 frames
        if self.seq_type == "P":
            return self.seq
        else:
            return ret_protein

    # instance method that returns GC content if seq_type is "N", else return 0
    def gc_content(self):
        gc_count = self.seq.count("G") + self.seq.count("C")
        if self.seq_type == "N":
            return float(gc_count / len(self.seq))
        return 0

    def __len__(self):
        return len(self.seq)

    # concatenates sequences of the same type, else creates empty sequence of type "N"
    def __add__(self, other):
        while isinstance(self, Sequence) and isinstance (other, Sequence):
            if self.seq_type == other.seq_type:
                return self.seq + other.seq

            else:
                return Sequence()

        return Sequence()


    # ensures both objects are Sequence instances before comparing them
    def __eq__(a, b):
        if (not isinstance(a, Sequence)) or (not isinstance(b, Sequence)):
            return False

        while isinstance(a, Sequence) and isinstance(b, Sequence):
            for i in a.translateSeq():

                if i in b.translateSeq():
                    return True

                return False

    def __str__(self):
        return self.seq

  # init method takes two optional parameters, seq and seq_type
    def __init__(self, seq ="", seq_type="N"):
        self.seq = seq
        self.seq_type = seq_type

        # if sequence is provided, but no type, checks sequence type with helper method
        if seq != "" and seq_type == "":
            self.seq_type = "P" if Sequence.sequence_type(seq) == False else "N"

        # if sequence and type are provided, and nucleotide type is specified, checks sequence
        # type with helper method. This allows user to override and specify a protein type for a
        # string containing only "AGTC"
        if seq != "" and seq_type != "":
            if seq_type == "N":
                self.seq_type = "P" if Sequence.sequence_type(seq) == False else "N"

# Creates a Sequence instance for each entry, and returns list of all Sequence instances.
def sequenceList():
    inst_list = list()
    while True:
        user = input("Enter a nucleotide sequence:")
        if user == "":
            break

        seq_inst = Sequence(user)
        inst_list.append(seq_inst)

    return inst_list


# takes list of Sequence instances, and an integer specifying kmer length. Returns dictionary
# with all possible kmers, where the kmer is the key and the number of occurences is the value.
def kmerDict(seq_instances, size):
    ret_dict = {}
    kmer_length = size - 1

    for inst in seq_instances:
        if inst.seq_type == "N": # only considers nucleotide strings

            for start_pos in range(len(inst.seq) - kmer_length):
                substring = inst.seq[start_pos:start_pos+size]
                ret_dict.setdefault(substring, 0)
                ret_dict[substring] += 1

    return ret_dict

seq_instances = sequenceList()


# test with all 3mers, 4mers, and 5mers
three_mer = kmerDict(seq_instances, 3)
print(three_mer)

four_mer = kmerDict(seq_instances, 4)
print(four_mer)

five_mer = kmerDict(seq_instances, 5)
print(five_mer)
