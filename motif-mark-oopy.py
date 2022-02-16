#!/usr/bin/python3

import cairo
import argparse
import itertools


class Motif:
    """
    Pass in a nucleotide motif that will be responsible for generating all of the possible
    combinations.
    Takes into account
    Y Pyrimidines
    R Purines
    """
    def __init__(self, mot_string):
        self.motif = mot_string.upper()
        self.purines = ['a', 'g']
        self.pyrimidines = ['c', 't']
        self.combinations = list()
        self.run = True
        while self.run == True:
            if 'Y' in self.motif
                loc = [i for i, ltr in enumerate(self.motif) if ltr == 'Y']
                for y_ent in loc:
                    for ele in self.pyrimidines:
                        self.combinations.append(self.motif[:y_ent] + ele + self.motif[y_ent+1:])
                if 'R' not in self.motif:
                    self.run = False
            if 'R' in self.motif:
                loc = [i for i, ltr in enumerate(self.motif) if ltr == 'R']
                for r_ent in loc:
                    for ele in self.purines:
                        self.combinations.append(self.motif[:r_ent] + ele + self.motif[r_ent+1:])
            else:
                self.combinations.append(self.motif)


class Gene:
    """
    Add in Exons and intron information
    """
    pass


class Cairo:
    """
    Take in the Motif and Gene object and then draw it.
    """
    pass
