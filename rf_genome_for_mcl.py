from Bio import SeqIO
import glob
import sys

genome_list = sys.argv[1]
species_list = open(genome_list, "r")

for species in species_list:
    species = species.strip()
    genome = species + ".fa"
    genome_records = SeqIO.parse(genome, "fasta")

    count = 1
    reformatted_description = {}
    for record in genome_records:
        description = record.description
        species = species.replace("_", ".")
        new_description = species + "_" + str(count)
        reformatted_description[new_description] = str(record.seq)
        count += 1 

    new_file_name = species + "_rf.fa"
    file = open(new_file_name, "w")
    for desc, seq in reformatted_description.items():
        line = ">" + desc + "\n" + seq + "\n"
        file.write(line)
    file.close()
