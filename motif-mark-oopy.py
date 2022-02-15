#!/usr/bin/python3

import cairo


class Motif:
    """
    Pass in a nucleotide motif that will be responsible for generating all of the possible
    combinations.
    Takes into account
    Y Pyrimidines
    R Purines
    """
    __init__(self, mot_string):
        self.motif = mot_string.upper()
        self.purines = ['a', 'g']
        self.pyrimidines = ['c', 't']
        self.combinations = list()
        if 'Y' in self.motif:
            loc = [i for i, ltr in enumerate(self.motif) if ltr == 'Y']
            for y_ent in len(loc):
                # Need to fix this part still
                self.combinations.append(self.motif[loc[y_ent] = 'None')
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
