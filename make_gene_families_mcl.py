import os
import sys
import glob
from Bio import SeqIO

def check_data(mcl):

    """ Function to check MCL data """

    print("\n*****************************************************************")
    print("\nChecking MCL data")

    number_families = 0
    f = open(mcl, "r")
    for l in f:
        number_families += 1

    print("Number of family clusters identified:", number_families)

def create_sequence_dictionary(data):

    """ Function to create sequence database from fasta file """

    print("\n*****************************************************************")
    print("\nCreating sequence database from fasta records")

    seq_dict = {}
    number_of_records = 0
    records = SeqIO.parse(data, "fasta")

    for seq_record in records:
        number_of_records += 1
        description = seq_record.description
        sequence = seq_record
        seq_dict[description] = sequence

    print("\nDatabase constructed")
    print("\nNumber of records in dataset:", number_of_records)

    return seq_dict

def create_gene_family_dictionary(mcl):

    """ Function to convert mcl output to dictionary """

    print("\n*****************************************************************")
    print("\nCreating gene family database from MCL output")

    count = 0
    mcl_records = 0
    gene_family_dict = {}
    f = open(mcl, "r")
    for l in f:
        count+=1
        gene_family = l.split()
        mcl_records += len(gene_family)
        gene_family_dict["Gene_family_"+str(count)] = gene_family

    print("\nGene family database constructed")
    print("\nNumber of gene families from MCL:",count)
    print("\nNumber of records clustered by MCL:",mcl_records)

    return gene_family_dict

def write_cluster_fasta_files(seq_dict, gene_family_dict):

    """ Function to write gene family clusters to fasta files """

    print("\n*****************************************************************")
    print("\nWriting fasta files for MCL clusters")

    number_of_gene_families = len(gene_family_dict)
    progress = 0
    pct_progress_been = set()
    for name, gene_family in gene_family_dict.items():
        progress +=1
        v_gene_f = {}
        for gene in gene_family:
            record = seq_dict[gene]
            v_gene_f[record.description] = str(record.seq)

        file_name = str(name) + ".fa"
        f = open(file_name,"w")
        for description, sequence in v_gene_f.items():
            line = (">"+ description +"\n" + sequence + "\n")
            f.write(line)
            pct_progress = int(round((progress/number_of_gene_families)*100,1))
            if pct_progress not in pct_progress_been:
                line = int(pct_progress/2) * ("*")
                print(line)
            pct_progress_been.add(pct_progress)
        f.close()

    print("\nAll MCL clusters have now been written to fasta files")

def back_checks(seq_dict, gene_family_dict):

    """ Function to back check all the analysis """

    print("\n*****************************************************************")
    print("\nPerforming back checks and identifying unclustered genes")

    written_records = []
    for fasta_file in glob.glob("*.fa"):
        gene_family = SeqIO.parse(fasta_file,'fasta')
        for record in gene_family:
            written_records.append(str(record.description))

    all_records = []
    for record, seq in seq_dict.items():
        all_records.append(record)

    unclustered_records =  list(set(all_records)-set(written_records))

    print("\nUnclustered records = {}".format(len(unclustered_records)))
    print("\nWriting unclusterd_records to text file")

    f = open("unclusterd_records.txt", "w")
    for seq in unclustered_records:
        f.write(seq +"\n")
    f.close()

    print("\nBack Checks finished")

    return(unclustered_records,all_records)

def clean():

    """ Function to clean directory """

    print("\n*****************************************************************")
    print("\nCleaning directory")

    os.mkdir("MCL_gene_families")

    total_gene_families = 0
    for fasta_file in glob.glob("*.fa"):
        total_gene_families += 1

    pct_progress_been = set()
    progress = 0
    for fasta_file in glob.glob("*.fa"):
        progress += 1
        command = "mv " + fasta_file + " MCL_gene_families"
        os.system(command)
        pct_progress = int(round((progress/total_gene_families)*100,1))
        if pct_progress not in pct_progress_been:
            line = ((int(pct_progress/2)) * ("*")) + " %" + str(pct_progress)
            print(line)
        pct_progress_been.add(pct_progress)

    print("\nAll gene familes have been moved ")

    os.mkdir("input_files")
    command = "mv " +  mcl + " input_files"
    os.system(command)

    command = "mv " + data + " input_files"
    os.system(command)

    os.mkdir("unclustered_records")
    command = "mv unclusterd_records.txt unclustered_records"
    os.system(command)

    print("\nAll cleaning done")
    print("\nGene families are in the MCL_gene_families directory")
    print("\nThere were {} gene families created".format(total_gene_families))
    print("\n*****************************************************************")

def main():
    check_data(mcl)
    seq_dict = create_sequence_dictionary(data)
    gene_family_dict = create_gene_family_dictionary(mcl)
    write_cluster_fasta_files(seq_dict,gene_family_dict)
    unclustered_records, all_records = back_checks(seq_dict,gene_family_dict)
    clean()

if __name__ == "__main__":
    mcl = sys.argv[1]
    data = sys.argv[2]
    main()
