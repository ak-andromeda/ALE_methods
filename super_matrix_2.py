import os
import sys
import glob
import time
import subprocess
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import matplotlib.pyplot as plt


def count_AA_occurence(sequence, amino_acids_to_track='ACDEFGHIKLMNPQRSTVWY-'):

    """ Function to count amino acide occurences """

    amino_acids_per_train = 0
    amino_acids = dict.fromkeys(amino_acids_to_track, 0)

    for char in sequence:
        if char in amino_acids:
            amino_acids_per_train += 1
            amino_acids[char] += 1

    AA_percentages = {k:(v*100.0)/amino_acids_per_train for k,
                   v in amino_acids.items()}
    return AA_percentages

def make_species_dict(species_file):

    """ Function to make dictionary of species present """

    print("*******************************************************************")
    print("\nSpecies sequence database being constructed")

    species_dict = defaultdict(list)
    with open(species_file) as f:
        for l in f:
            line = l.rstrip()
            species = line.replace(" ", ".")
            species_dict[species] = []

    # Read in the alignments
    for alingment in glob.glob("*.aln"):
        inhandle = SeqIO.parse(alingment, "fasta")

        present_species = set()
        all_species = set()

        for record in inhandle:

            # Reformat description to be the same as search key
            description = str(record.description)
            sequence = str(record.seq)
            seq_length = len(sequence)
            description_search = description.replace(" ", "_")

            # Search description to see if species if present
            for key, val in species_dict.items():
                key_search = key
                all_species.add(key_search)
                if key_search in description_search:
                    species_dict[key_search].append(sequence)
                    present_species.add(key_search)

        # Fill in gaps for all species with no sequences.
        not_present = all_species - present_species
        for key_search in not_present:
            species_dict[key_search].append("-" * seq_length)

    print("\nSpecies sequence database constructed")

    return species_dict

def make_concatentation(species_dict):

    """ Function to make the concatenation """

    print("\n*****************************************************************")
    print("\nConcatentation being constructed")


    Concatenated_seqs = {}
    for key,val in species_dict.items():
        concat = "".join(val)
        Concatenated_seqs[key] = concat

    timestr = time.strftime("%d%m%y-%H%M")
    file_name = "Concatenate_" + timestr + ".fasta"
    ofile = open(file_name, "w")
    for key, val in Concatenated_seqs.items():
        ofile.write(">" + key + "\n" + val + "\n")
    ofile.close()

    print("\nConcatenation made and saved to:", file_name)

    return Concatenated_seqs

def analyze_gap_content(Concatenated_seqs):

    """ Function to analyse the alignment """

    print("\n*****************************************************************")
    print("\nConcatentation gap content being analysed")

    Species = []
    Sequences = []
    for key, val in Concatenated_seqs.items():
        Species.append(key)
        Sequences.append(val)
    data = {'Species': Species, 'Sequences': Sequences}
    df_seq = pd.DataFrame.from_dict(data)

    Species_AA_content = {}
    for species, seq in Concatenated_seqs.items():
        AA_content = count_AA_occurence(seq)
        Species_AA_content[species] = AA_content

    df_AA = pd.DataFrame.from_dict(Species_AA_content)
    df_AA = df_AA.transpose()
    df_AA = df_AA.reset_index()
    df_AA.rename(columns={ df_AA.columns[0]: "Species" }, inplace = True)

    df_concat = df_seq.merge(df_AA, on = "Species")
    df_concat.to_csv("df_concat.csv")

    print("\nResults of gap and amino acid content saved to df_concat.csv")
    return df_concat

def plot_gap_content(df_concat):

    """ Function to plot the results of gap analysis """

    print("\n*****************************************************************")
    print("\nPlots of gap content being produced")

    df_concat = df_concat.sort_values(by=['-'])
    sns.set(rc={'figure.figsize':(15,10)})
    sns.set_style("whitegrid", {"axes.grid":False})
    Gap_content = sns.barplot(x="Species", y="-", data=df_concat,
                              palette = "flare")
    Gap_content.set_title("Concatenation Gap Content")
    Gap_content.set_xlabel("Species")
    Gap_content.set_ylabel("% of gaps in alignment")
    plt.xticks(rotation=90)
    plt.tight_layout()
    figure = Gap_content.get_figure()
    figure.savefig('Concatenation_gap_content.png', dpi=400)

    print("\nPlots saved to: Concatenation_gap_content.png")

def make_conservative_concatentate(df_concat,Concatenated_seqs):

    """ Function to make concatentate without high gap species """

    print("\n*****************************************************************")
    print("\nMaking conservative concatenate - >50% gap species removed")

    timestr = time.strftime("%d%m%y-%H%M")

    df_conservative = df_concat[df_concat["-"] <= 50]
    conservative_species = list(df_conservative.Species)

    file_name = "Conservative_concatenate_" + timestr + ".fasta"
    ofile = open(file_name, "w")
    for species in conservative_species:
        ofile.write(">" + species + "\n" + Concatenated_seqs[species] + "\n")
    ofile.close()

    file_name = "Conservative_species_list_" + timestr + ".txt"
    ofile = open(file_name, "w")
    for species in conservative_species:
        ofile.write(species + "\n")
    ofile.close()

    problem_species = list(df_concat[df_concat["-"] >= 50].Species)
    file_name = "Species_removed_from_concat_" + timestr + ".txt"
    ofile = open(file_name, "w")
    for species in problem_species:
        ofile.write(species + "\n")
    ofile.close()

    print("\nConservative concatenate made species removed =", problem_species)
    return problem_species

def clean_up():

    """ Function to clean up directory """

    print("\n*****************************************************************")
    print("\nDirectory being cleaned")

    # Sort files
    os.mkdir("Results")
    os.mkdir("Concatenations")
    os.mkdir("Original_trimmed_orthogroups")
    os.mkdir("Super_matrix_py_script")

    my_commands = ["mv *.csv Results",
                   "mv *.txt Results"
                   "mv *.png Results",
                   "mv *.fasta Concatenations",
                   "mv *.aln Original_trimmed_orthogroups",
                   "mv *.py Super_matrix_py_script"]

    for my_command in my_commands:
        os.system(my_command)

def output(df_concat):

    """ Function to print a summary output """

    print("\n*****************************************************************")
    print("\nSummary")

    sequence = df_concat.at[0,'Sequences']
    print("Length of alignment =",len(sequence))
    print("Number of species = ", len(df_concat))
    print("\nResults are now in the the 'Results' directory")
    print("Both the concatentations are in the 'Concatenations'directory")
    print("The individual orthogroups are in 'Original_trimmed_orthogroups'")
    print("\nThank you for using super_matrix_2.py")

def main():
    species_dict = make_species_dict(species_file)
    Concatenated_seqs = make_concatentation(species_dict)
    df_concat = analyze_gap_content(Concatenated_seqs)
    problem_species = make_conservative_concatentate(df_concat,Concatenated_seqs)
    plot_gap_content(df_concat)
    clean_up()
    output(df_concat)

if __name__ == "__main__":
    species_file = sys.argv[1]
    main()
