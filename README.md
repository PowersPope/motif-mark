# Motif Mark using Cairo

## Introduction

Motif Mark is a tool that will visually represent motifs on a gene (introns/exons). We
utilize the power of PyCairo and create 1:1 realistic representations of the gene and 
location of the motifs. Written in Python 3.9.


## Set Up

You will need PyCairo installed to make this work. We have made it easy to install the 
dependencies if you have conda installed on your computer already. Just follow the command
below if you are within the repository directory directly. Then you are good to go!

```
conda env create -f env_cairo.yml
```


## Parameters

There are two parameters/files you that have to supply to the tool in order for it to work.
A fasta file with your intron/exons, one exon per fasta line. Also, you will need to supply 
a motif file. Which can be a text file, however, what is necessary is that each motif is on
a newline. Like the example below.

Motif Example:

```
ygcy
GCAUG
catag
YYYYYYYYYYY
```

When calling the python file you can do it as below:

```
python motif-mark-oopy.py -f <fasta file> -m <motif file>
```

You can see that you supply the fasta file with `-f` and the motif file with `-m` .


This will output a .png file that has the same name as the file passed in the `-f` parameter.


