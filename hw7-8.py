#! /usr/bin/env python3

import os

all_species = dict()

try:
    with open("hw7-input.sam", "r") as file_handle:
        for line in file_handle:
            line = line.rstrip()

            # Skip SAM header
            if line.startswith("@"):
                continue

            # Extract fields
            fields = line.split("\t")

            # Skip unmapped
            if fields[1] == '4':
                continue

            # Parse and record species
            species = fields[2].split("|")[2].split("_")[1]
            all_species.setdefault(species, 0)
            all_species[species] += 1
except IOError as test:
    print("IO Error!",test)

# Report results
for species in all_species:
    print("{0}: {1}".format(species, all_species[species]))
