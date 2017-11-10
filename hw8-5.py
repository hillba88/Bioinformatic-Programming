#! /usr/bin/env python3

import Bio.SeqIO
import Bio.Seq
import argparse

# create instance of argparse, with four positional arguments
parser = argparse.ArgumentParser()
parser.add_argument("path", help="input file path")
parser.add_argument("out", help="output file path")
parser.add_argument("frame", type=int, help="reading frame")
parser.add_argument("table", help="translation table")

args = parser.parse_args()


try:
    records = []
    # iterate through input file and translate each sequence
    for record in Bio.SeqIO.parse(args.path, "fasta"):
        strand = record.seq if args.frame < 3 else record.seq.reverse_complement()

        trans_seq = ""
        # translates forward or reverse complement, depending on frame and table arguments
        for pos in range(args.frame % 3, len(strand), 3):
            codon = strand[pos:pos+3]
            trans_seq += codon.translate(table=args.table)

        # creates new SeqRecord for each translation and appends to a list
        cur_record = Bio.SeqRecord.SeqRecord(id=record.id, seq=record.seq, description='')
        records.append(cur_record)

    # writes list of translated SeqRecords to output file in fasta format
    Bio.SeqIO.write(records, args.out, "fasta")

except IOError:
    print("IO Error!")
