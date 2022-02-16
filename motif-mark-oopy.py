#!/usr/bin/python3

import cairo
import argparse
import itertools as it
import re


def get_args():
    args = argparse.ArgumentParser(description="Pass in Gene Fasta and Motif files")
    args.add_argument(
        "-f", "--file", help="Path to Gene Fasta file", required=True, type=str
    )
    args.add_argument(
        "-m", "--motif", help="Path to Motif file", required=True, type=str
    )
    return args.parse_args()


args = get_args()


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

    def search_gene(self, pre, exon, post):
        """Search through a gene object and find all of the motif matches"""
        pre_dict = dict()
        exon_dict = dict()
        post_dict = dict()
        for ent in self.combinations:
            pre_matches = re.finditer(ent, pre)
            exon_matches = re.finditer(ent, exon)
            post_matches = re.finditer(ent, post)
            if [stu for stu in pre_matches] != list():
                pre_dict[ent] = pre_matches
            elif [stu for stu in exon_matches] != list():
                exon_dict[ent] = exon_matches
            elif [stu for stu in post_matches] != list():
                post_dict[ent] = post_matches
        return (pre_dict, exon_dict, post_dict)


class Gene:
    """
    Pass in the lines that consitute a fasta entry from a fasta file. Grab the Introns (lowercase)
    and Exons (uppercase) and put them into separate fields of pre Exon intron and post
    Exon intron.
    """

    def __init__(self, fasta_entry_lines, motif_dict):
        self.header = fasta_entry_lines[0]
        self.gene = fasta_entry_lines[1]
        self.parts = re.findall(r"([a-z]+)([A-Z]+)([a-z]+)", self.gene)
        self.pre_exon = self.parts[0][0]
        self.exon = self.parts[0][1]
        self.post_exon = self.parts[0][2]
        self.gene_len = len(self.gene)
        self.motif_matches = dict()

        for motif in motif_dict:
            self.motif_matches[motif] = motif_dict[motif].search_gene(
                self.pre_exon, self.exon, self.post_exon
            )

    def show_matches(self):
        return self.motif_matches


class Cairo:
    """
    Take in the Motif and Gene object and then draw it.
    """

    pass


####### Logic

# init dictionary full of all motifs
motif_dict = dict()

# load in motif file
with open(args.motif, "r") as motif_file:
    for line in motif_file:
        # Strip the new line characters off of it
        clean = line.strip("\n")
        # Add each line as the key and then create a motif object from the line
        motif_dict[clean] = Motif(clean)


print(motif_dict)

# Init global variables that will be used in the for and with loop below
num = 0
gene_dict = dict()
entry_holder = list()
unwrapped_line = ""

# Open the gene file
with open(args.file, "r") as file:
    for line in file:
        # Check to see that it is a header line
        if r">" in line:
            # Just add the line if unwrapped_line == ''
            if unwrapped_line == "":
                entry_holder.append(line.strip("\n").lstrip(">"))
            else:
                # Add entry to dict
                num += 1
                entry_holder.append(unwrapped_line)
                print("Making Object:", entry_holder)
                gene_dict[num] = Gene(entry_holder, motif_dict)
                # Reset the holding variables
                unwrapped_line = ""
                entry_holder.clear()
                entry_holder.append(line.strip("\n").lstrip(">"))
        else:
            unwrapped_line += line.strip("\n")
    # Add the last entry to the dictionary
    num += 1
    entry_holder.append(unwrapped_line)
    gene_dict[num] = Gene(entry_holder, motif_dict)

print(gene_dict)


print(gene_dict[1].show_matches())

