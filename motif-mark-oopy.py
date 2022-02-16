#!/usr/bin/python3

import cairo
import argparse
import itertools as it
import re


class Motif:
    """
    Pass in a nucleotide motif that will be responsible for generating all of the possible
    combinations.
    Takes into account
    Y Pyrimidines
    R Purines
    Returns a list of all possible combinations
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
            self.combinations.append("".join([i for i in combos]))

    def combos(self):
        """Return all possible motif combinations in a list"""
        return self.combinations

    def search_gene(self, gene):
        """Search through a gene object and find all of the motif matches"""
        pre = gene.pre_exon
        exon = gene.exon
        post = gene.post_exon
        pre_dict = dict()
        exon_dict = dict()
        post_dict = dict()
        for ent in self.combinations:
            pre_matches = re.finditer(ent, pre)
            exon_matches = re.finditer(ent, exon)
            post_matches = re.finditer(ent, post)
            if [stu for stu in pre_matches] != list():
                pre_dict[ent] = pre_matches
            else:
                pre_dict[ent] = "None"
            if [stu for stu in exon_matches] != list():
                exon_dict[ent] = exon_matches
            else:
                exon_dict[ent] = "None"
            if [stu for stu in post_matches] != list():
                post_dict[ent] = post_matches
            else:
                post_dict[ent] = "None"
        return pre_dict, exon_dict, post_dict


class Gene:
    """
    Pass in the lines that consitute a fasta entry from a fasta file. Grab the Introns (lowercase)
    and Exons (uppercase) and put them into separate fields of pre Exon intron and post
    Exon intron.
    """

    def __init__(self, fasta_entry_lines):
        self.header = fasta_entry_lines[0]
        self.gene = fasta_entry_lines[1]
        self.parts = re.findall(r"([a-z]+)([A-Z]+)([a-z]+)", self.gene)
        self.pre_exon = self.parts[0][0]
        self.exon = self.parts[0][1]
        self.post_exon = self.parts[0][2]
        self.gene_len = len(self.gene)


class Cairo:
    """
    Take in the Motif and Gene object and then draw it.
    """

    pass


##########

mot = Motif("ygtcrcty")
gene = Gene(["test", "cgtcgcttctgattatgGTCATAGTCCATATgtacgtcgctttctagcgtcgctt"])

print(mot.combos())
print(gene.pre_exon)
print(gene.exon)
print(gene.post_exon)
for amp in mot.search_gene(gene):
    print(amp)
