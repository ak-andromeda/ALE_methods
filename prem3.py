#!/usr/bin/python3
import os
import sys
import glob
from Bio import SeqIO
from Bio.Seq import Seq

def make_species_set(species_file):

    """ Function to make a species database """

    species_set = set()
    with open(species_file) as f:
            for l in f:
                species = l.rstrip()
                species_set.add(species)

    print("***************************************\n")
    print("Total number of species under analysis: ", len(species_set))
    return species_set

def calculate_total_orthogroups():

    """ Function to calculate total number of orthogroups """

    total_orthogroups = 0
    for orthogroup in glob.glob("*.fa"):
        total_orthogroups += 1

    return total_orthogroups

def make_single_copy_families(total_orthogroups, species_set):

    """ Function to remove paralogs and replace with gaps """

    orthogroup_count = 0
    species_rep = {}

    for orthogroup in glob.glob("*.fa"):
        orthogroup_count += 1
        print("Orthogroup under analysis: ", orthogroup)
        print("Orthologs analysed: ", orthogroup_count, "/", total_orthogroups)
        inhandle = SeqIO.parse(orthogroup, "fasta")

        species_count = {}
        for species in species_set:
            species_count[species] = 0

        all_descriptions = {}
        for record in inhandle:
            record_description = record.description.replace(" ", "_")
            all_descriptions[record_description] = str(record.seq)
            for species in species_set:
                if species in record_description:
                    species_count[species] += 1

        sc_orthgroup_species = set()
        mc_orthgroup_species = set()
        nc_orthgroup_species = set()
        removed_species = set()

        for key, val in species_count.items():
            if val == 1:
                sc_orthgroup_species.add(key)
            elif val > 1:
                mc_orthgroup_species.add(key)
                removed_species.add(key)
            else:
                nc_orthgroup_species.add(key)
                removed_species.add(key)

        ofile = open(orthogroup[:9] + "_sc.faa", "w")
        for description, seq in all_descriptions.items():
            for species in sc_orthgroup_species:
                if species in description:
                    ofile.write(">" + description + "\n" + seq + "\n")
        ofile.close()

        Total_species = len(species_set)
        Species_with_sc = len(sc_orthgroup_species)
        Species_with_mc = len(mc_orthgroup_species)
        Species_with_nc = len(nc_orthgroup_species)

        sequences = 0
        sc_orthogroup = SeqIO.parse(orthogroup[:9] +"_sc.faa", "fasta")
        for sequence in sc_orthogroup:
            sequences += 1
        if sequences == Species_with_sc:
            print("sc orthofile correctly printed")
        else:
            print("Error with code")

        species_representation = (Species_with_sc/Total_species)*100
        o_file = orthogroup[:9] + "_sc.faa"
        species_rep[o_file] = species_representation

    return species_rep

def move_orthogroups(total_orthogroups,species_rep, cutoff):

    """ Function to move files based on % species representation """

    more_dir = "Orthogroups_more_than_%" + str(cutoff)
    less_dir = "Orthogroups_less_than_%" + str(cutoff)

    os.mkdir(more_dir)
    os.mkdir(less_dir)

    ofile = open("species_rep.csv","w")
    more_files = 0
    completed = 0
    for key, val in species_rep.items():
        ofile.write(key + "," + str(val) + "\n")
        completed += 1
        if val >= int(cutoff):
            more_files += 1
            my_command =  "mv " + key + " " + more_dir
            os.system(my_command)
            line1 = "\nCurrent single copy families:{}".format(more_files)
            line2 = "Searched: %" + str((completed/total_orthogroups)*100)
            print(line1,line2)

        if val <= int(cutoff):
            my_command =  "mv " + key + " " + less_dir
            os.system(my_command)

    return(more_files)

def clean_up():

    """ Move original files """

    print("\n*****************************************************************")
    print("Cleaning up - all files being moved to labled sub-directory")
    os.mkdir("Original_orthogroups")
    os.mkdir("Results")

    my_c = """for file in *.fa; do \n mv "$file" Original_orthogroups \n done"""
    os.system(my_c)

    my_c = "mv *.csv Results"
    os.system(my_c)

def output(more_files,cutoff):

    print("\n*****************************************************************")
    print("\nprem3 has finished")
    line1 = "All single copy families above %{}".format(cutoff)
    line2 = "have been moved to their own directory"
    print(line1,line2)
    print("{} single copy files have been generated".format(more_files))
    print("Thanks for using prem3")
    print("\n*****************************************************************")

def main():
    species_set = make_species_set(species_file)
    total_orthogroups = calculate_total_orthogroups()
    species_rep = make_single_copy_families(total_orthogroups,species_set)
    more_files = move_orthogroups(total_orthogroups,species_rep, cutoff)
    clean_up()
    output(more_files,cutoff)

if __name__ == "__main__":
    species_file = sys.argv[1]
    cutoff = sys.argv[2]
    main()
