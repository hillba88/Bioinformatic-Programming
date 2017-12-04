#! /usr/bin/env python3

import random 

class Lineage:
    def __init__(self, length):
        self.length = random.randit(int(0.75*length), int(1.25*length))
        self.genomes = []
        self.genomes.append(gen_seq(self.length))
        self.active = 0

    def divide(self, prob):
        new_seq = ""
        for nuc in self.genomes[-1]:
            if random.random() < prob:
                new_seq += random.choice("ACTG")
            else:
                new_seq += nuc
    
    def set_active(self, active):
        if active >= len(self.genomes):
            active = len(self.genomes) - 1

        self.active = active

    def __str__(self):
        return self.genomes[self.active]


