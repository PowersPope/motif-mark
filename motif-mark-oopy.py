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
        if "u" in self.motif:
            self.motif = self.motif.replace("u", "t")
        self.purines = "[ag]"
        self.pyrimidines = "[ct]"
        self.motif_regex = mot_string.lower()

        if "y" in self.motif:
            self.motif_regex = self.motif.replace("y", self.pyrimidines)
        elif "r" in self.motif_regex:
            self.motif_regex = self.motif_regex.replace("r", self.purines)

    def combos(self):
        """Return all possible motif combinations in a list"""
        return self.motif_regex

    def search_gene(self, pre, exon, post):
        """Search through a gene object and find all of the motif matches"""
        # init holding dicts
        match_dict = dict()
        # grab the length of each item
        pre_len = len(pre)
        exon_len = len(exon)
        post_len = len(post)
        preAndexon = exon_len + pre_len
        # Find the particular items of interest
        pre_matches = re.finditer(self.motif_regex, pre)
        exon_matches = re.finditer(self.motif_regex, exon)
        post_matches = re.finditer(self.motif_regex, post)

        for item in pre_matches:
            if item.group() in match_dict:
                match_dict[item.group()].append(item.span())
            else:
                match_dict[item.group()] = [item.span()]

        for item in exon_matches:
            new_span = (item.span()[0] + pre_len, item.span()[1] + pre_len)
            if item.group() in match_dict:
                match_dict[item.group()].append(new_span)
            else:
                match_dict[item.group()] = [new_span]

        for item in post_matches:
            new_span = (item.span()[0] + preAndexon, item.span()[1] + preAndexon)
            if item.group() in match_dict:
                match_dict[item.group()].append(new_span)
            else:
                match_dict[item.group()] = [new_span]

        return match_dict


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

    def __init__(self, gene_dict):
        self.width = max([gene_dict[n].gene_len for n in gene_dict]) + 200
        self.height = 600


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


# print(motif_dict)

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

# print(gene_dict)
for num in range(1, len(gene_dict) + 1):
    print(gene_dict[num].show_matches())

# for k in gene_dict[1].show_matches():
# for i in gene_dict[1].show_matches()[k]:
# print(i)

