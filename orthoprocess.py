# Scrip to remove specific species from alignments.
# Written by Brogan Harris - 15/08/2019

# Imports and Libraries required #
import os
import sys
import glob
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq


os.mkdir("four_plus_seq")
os.mkdir("minus_four_seq")

count = 0
for orthogroup in glob.glob("*fa*"):
    count += 1

analysed_count = 0 

# Read in the alignments
for orthogroup in glob.glob("*fa*"):
    ortho_seqs = SeqIO.parse(orthogroup, "fasta")

    num_seqs = 0 
    for record in ortho_seqs:
        num_seqs += 1
    
    if num_seqs >= 4:
        my_command = "mv " + orthogroup + " four_plus_seq"
        output  = subprocess.check_output(my_command, shell=True,
                                  stderr=subprocess.STDOUT)
        print(orthogroup)
    
    else: 
        my_command = "mv " + orthogroup + " minus_four_seq"
        output  = subprocess.check_output(my_command, shell=True,
                                  stderr=subprocess.STDOUT)
        print(orthogroup)
    
    analysed_count += 1

    print("orthogroups analysed = ", analysed_count, "/", count)



