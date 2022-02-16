#!/usr/bin/python3

import cairo
import argparse
import itertools as it


class Motif:
    """
    Pass in a nucleotide motif that will be responsible for generating all of the possible
    combinations.
    Takes into account
    Y Pyrimidines
    R Purines
    """
    def __init__(self, mot_string):
        self.motif = mot_string.lower()
        self.purines = ["a", "g"]
        self.pyrimidines = ["c", "t"]
        self.combinations = list()
        self.dict_comb = dict()
        for i, j in enumerate(self.motif):
            if j == "y":
                self.dict_comb[i] = self.pyrimidines
            elif j == "r":
                self.dict_comb[i] = self.purines
            else:
                self.dict_comb[i] = [j]
        for combos in it.product(*(self.dict_comb[k] for k in self.dict_comb)):
            self.combinations.append(''.join([i for i in combos]))
    def show_combos(self):
        return self.combinations


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
