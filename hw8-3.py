#! /usr/bin/env python3
import csv
import subprocess
import argparse

# create instance of argparse, add two argument, one for input file and
# one for output file
parser = argparse.ArgumentParser()
parser.add_argumnet("--reads", help="path for forward reads")
parser.add_argument("--output", help="ouput basename")

args = parser.parse_args()

# run two paladin prepare and align, using input and output arguments
paladin_prepare = subprocess.run("paladin prepare -r1", shell=True,
                                 stdout=subprocess.PIPE)
paladin_align = subprocess.run("paladin align uniprot_sprot.fasta.gz {0} -t 10 -o {1}".format(args.reads, args.output),
                               shell=True, stdout=subprocess.PIPE)


try:
# reads SAM file as tab delimited file
    with open("{0}.sam".format(args.output), "r", newline="") as my_file:
        reader = csv.reader(my_file, delimiter = "\t")
        mapped_reads = []
        # iterates through each line and conditionally appends mapping quality to a list
        for line in reader:
            if line[1] != 4 and "@" not in line[0]: # only appends if reads are mapped, skips headers
                mapped_reads.append(line)

        # loop through mapped reads and append reference sequence names to a list
        ref_names = []
        for item in mapped_reads:
            if item[2] != "*":
                ref_names.append(item[2])

        # loop through reference sequence names, split, and append KBID to a list
        kbid = []
        for item in ref_names:
            item = item.split("|")
            kbid.append(item[2])

        # loop through KBIDs and add species IDs to a list
        species = []
        for item in kbid:
            item = item.split("_")
            species.append(item[1])

        # populate a dictionary with the species ID and how many times it appears in the SAM file
        species_dict = {}
        for item in species:
            species_dict.setdefault(item, 0)
            species_dict[item] += 1

        # print species ID and count, one per line
        for key,val in species_dict.items():
            print(key,":",val)

except csv.Error:
    print("Invalid file path")
