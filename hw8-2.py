#! /usr/bin/env python3

import subprocess
import argparse

# # create instance of argparse, add two argument, one for input file and
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
