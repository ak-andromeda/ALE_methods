# Standard library packages
import io
import os
import sys
from collections import Counter

# Import Biopython modules to interact with KEGG
from Bio import SeqIO
from Bio.KEGG import REST
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas


def read_in_kegg_data(kegg_data,node):

    kegg_dictionary = {}
    file = open(kegg_data, "r")

    for l in file:
        l = l.strip()
        l = l.split("\t")

        if len(l) >= 2:
            gene_family = l[0]
            kegg = l[1]
            kegg_dictionary[gene_family] = kegg

    file_name = node + "_kegg_dict.csv"
    file = open(file_name, "w")
    for key, val in kegg_dictionary.items():
        line = key + "," + str(val) + "\n"
        file.write(line)
    file.close()

    return kegg_dictionary


def find_kegg_pathway(kegg_dictionary,node):

    pathways = []
    pathway_dictionary = {}
    for gene_family, kegg in kegg_dictionary.items():
        result = REST.kegg_get(kegg)
        for l in result:
            if l.startswith("P"):
                pathway = l.split("  ")[3].strip()
                pathway = pathway.replace(",", "")
                print(pathway)
                pathway_dictionary[gene_family] = pathway
                pathways.append(pathway)

    d = Counter(pathways)
    print(d)

    file_name = node + "_pathway_freq.csv"
    file = open(file_name, "w")
    for key, val in d.items():
        line = key + "," + str(val) + "\n"
        file.write(line)
    file.close()


def convert_to_df(node):

    keggs = []
    frequency = []

    file_name = node + "_pathway_freq.csv"
    file = open(file_name,"r")
    for line in file:
        line = line.strip().split(",")
        keggs.append(line[0])
        frequency.append(line[1])

    nodes = [node] * len(keggs)
    data = {"Node":nodes,"Kegg_pathway":keggs,"Frequency":frequency}
    df = pd.DataFrame.from_dict(data)

    df_name = node + "_pathway_freq_df.csv"
    df.to_csv(df_name)

    return df

def main():
    nodes = ["23"]
    for node in nodes:
        kegg_data = node + "_user_ko.txt"
        kegg_dictionary = read_in_kegg_data(kegg_data,node)
        find_kegg_pathway(kegg_dictionary,node)
        df = convert_to_df(node)



if __name__ == "__main__":
    main()
