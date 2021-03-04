import sys
import glob
import os
import time
import subprocess
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


def usage():

    """ Function to print usage """

    print ("""
############################################################################
Usage:

This script requires the following libraries to be installed:
sys, glob, os, time, subprocess, pandas, numpy, seaborn, biopython
and matplotlib. Version >=0.9 for seaborn is required for chart palettes. 
If you have an error regarding the use of the 'flare' palette, remove this 
argument from the plot_results() function. 

The script requires two text files named "roots_to_test.txt", and
"species_list.txt". The roots_to_test.txt file should
contain the names of the directories of the rec_files.
The species_list.txt file should contain all the species in the
dataset and should match the formatting of species tree used
in the reconciliation analysis.

The script also requires two command line arguments:

1) The maximum-likelihood root determined by an AU test
2) The metric you want to test: LS, DS or TS.

High species representation is currently set to %50 - this can be manually
changed by altering the make_high_rep_df function.

Run the script in the directory where all the directories containing 
the uml_rec files for each candidate root are location

Approximate time for a dataset of 10,000 uml_rec = 10 minutes.

############################################################################
    """)


class cd:

    """ Context manager for changing the current working directory """

    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

def get_count(ML_root):

    """ Function to find number of REC files for each root """

    new_path = os.getcwd() + "/" + ML_root
    with cd(new_path):
        command = "ls *uml_rec | wc -l"
        sp_output = subprocess.check_output(command, shell=True,
                                            stderr=subprocess.STDOUT)
        count = sp_output.strip()
        count = int(count)

    print("\nGene family rank is now being calculated")
    print("\n*******************************************","\n")
    print("Number of uml_rec files in each root:", count, "\n")
    print("\n*******************************************","\n")

    return count

def get_species_rep(ML_root):

    """ Function to calculate species representation """

    print("Species representation for each uml_rec file being calculated")
    file = open("species_list.txt", "r")
    species_list = []
    for i in file:
        species = i.strip()
        species_list.append(species)

    rtrb = ML_root
    new_path = os.getcwd() + "/" + rtrb
    species_rep_dictionary = {}
    with cd(new_path):
        current_dir = os.getcwd()
        for file_name in glob.glob("*.uml_rec"):
            species_count = 0
            rec_file = open(file_name, "r")
            for i, l in enumerate(rec_file):
                if i == 13:
                    species_tree = l
                    for species in species_list:
                        if species in species_tree:
                            species_count +=1
            species_rep = (species_count/len(species_list))*100
            species_rep_dictionary[file_name] = species_rep
            rec_file.close()

    timestr = time.strftime("%d%m%y")
    file_name = "species_repsentation_" + timestr + ".txt"
    file = open(file_name, "w")
    for key, val in species_rep_dictionary.items():
        line = key + "," + str(val) + "\n"
        file.write(line)
    file.close()

    print("\nSpecies representation calculated and saved in the text file:",
           file_name)

    print("\n*******************************************","\n")

    return(species_rep_dictionary)

def get_DTL(rec_file):

    """ Function to extract the loss/speciation ratio from rec_file """

    ratio_dict = {}
    for l in rec_file:
        if l.startswith("Total"):
            results = l.strip()
            result_list = results.split(" ")
            DTLS = result_list[1]
            DTLS_list = DTLS.split("\t")[1:]

            d = float(DTLS_list[0])
            t =  float(DTLS_list[1])
            l = float(DTLS_list[2])
            s = float(DTLS_list[3])

            if d == 0:
                d = 1

            if t == 0:
                t = 1

            if l == 0:
                l = 1

            if s == 0:
                s = 1


            ds = d/s
            ts = t/s
            ls = l/s

            ratio_list = [ds,ts,ls]

    return ratio_list

def sort_dict(dict):

    """ Function to sort dictionary for lowest to highest """

    sorted_dict = {k:v for k,v in sorted(dict.items(),key=lambda item: item[1])}
    return(sorted_dict)

def make_sorted_metric_dict(rank_by, DTLS_dict):

    """Function make a specific dictionary from the DTLS dictonary"""

    ranks = {"DS":0,"TS":1,"LS":2}

    rank_dict = {}
    if rank_by == "DS":
        index = ranks[rank_by]
        for key, val in DTLS_dict.items():
            rank_dict[key] = val[index]

    elif rank_by == "TS":
        index = ranks[rank_by]
        for key, val in DTLS_dict.items():
            rank_dict[key] = val[index]

    elif rank_by == "LS":
        index = ranks[rank_by]
        for key, val in DTLS_dict.items():
            rank_dict[key] = val[index]

    else:
        print("Not a suitable ranking metric")

    sorted_values = sort_dict(rank_dict)

    return sorted_values

def get_rec_file_rank(sorted_values, count):

    """ Function to assign rank to orthogroup based on metric"""

    number_rec_files = int(count)
    rec_file_rank = {}
    for file_name, metric in sorted_values.items():
        rec_file_rank[file_name] = number_rec_files
        number_rec_files -=1

    return rec_file_rank

def get_log_likelihood(rec_file):

    """ Function to extract the loss/speciation ratio from rec_file """

    for l in rec_file:
        if l.startswith(">logl:"):
            results = l.strip()
            results = float(results[7:])

    return results

def run_DTL(count,ML_root,rank_by,species_rep):

    """ function to execute DTL calculation functions """


    print("DTL for each gene family are being calculated")

    rtrb = ML_root
    rtt  = open("roots_to_test.txt", "r")
    new_path = os.getcwd() + "/" + rtrb

    DTLS_dict = {}
    with cd(new_path):
        current_dir = os.getcwd()
        for file_name in glob.glob("*.uml_rec"):
            rec_file = open(file_name, "r")
            DTL_ratio_list = get_DTL(rec_file)
            DTLS_dict[file_name] = DTL_ratio_list

    sorted_values = make_sorted_metric_dict(rank_by,DTLS_dict)
    rec_file_rank = get_rec_file_rank(sorted_values,count)

    d = {'Root': [], 'Gene_family': [], "log_likelihood": [], rank_by: [] }
    df = pd.DataFrame(data=d)
    colums = list(df)

    for root in rtt:
        root = root.strip()
        print("\nRoot currently under analysis:", root)
        new_path = os.getcwd() + "/" + root
        with cd(new_path):
            current_dir = os.getcwd()
            print(current_dir)
            for file_name in glob.glob("*.uml_rec"):
                rec_file = open(file_name, "r")
                log_likelihood = get_log_likelihood(rec_file)
                species_present = species_rep[file_name]
                rank = rec_file_rank[file_name]
                data = {"Root":root,
                        "Gene_family":file_name,
                        "species_representation": species_present,
                        "log_likelihood":log_likelihood,
                        rank_by:rank}
                df = df.append(data, True)
                rec_file.close()

    return(df)

def write_results(df, rank_by):

    """ Function to write pandas to csv """

    timestr = time.strftime("%d%m%y")
    file_name = rank_by + "_" + timestr + "_ranked.csv"
    df.to_csv(file_name)
    print("\nResults have been written to",  file_name)
    print("\nSlice of dataframe shown below:")
    print(df.head())

def make_high_rep_df(df,rankby):

    """ make high species representation dataframe """

    timestr = time.strftime("%d%m%y")
    file_name = rank_by + "_" + timestr + "_high_rep_ranked.csv"

    df = df[df["species_representation"] >= 50]
    df.to_csv(file_name)

    print("\nResults from gene families with >=50% species representation" + \
          " have been written to:", file_name)

    print("\nSlice of dataframe shown below:")
    print(df.head())

def summed_likelihood(rank_by,count):

    """ Function to sum likelihoods based on rank for each root """

    print("\n*******************************************","\n")
    print("\nSummed likelihoods and ranks now being calculated" + \
          " for uml_rec files\n")

    LOG_EVERY_N = 100
    number_of_roots = 0
    rtt  = open("roots_to_test.txt", "r")
    for root in rtt:
        number_of_roots += 1
    print("Number of roots to sum likelihood for:", number_of_roots)
    rtt.close()

    timestr = time.strftime("%d%m%y")
    file_name = rank_by + "_" + timestr + "_ranked.csv"
    df = pd.read_csv(file_name)
    d = {rank_by: [], 'Root': [], "summed_likelihood": [] }
    df_totals = pd.DataFrame(data=d)

    rtt  = open("roots_to_test.txt", "r")
    progress = 0
    for root in rtt:
        for rank in range(0,int(count),1):
            root = root.strip()
            df_r = df[df["Root"] == root]
            df_l = df_r[df_r[rank_by]>= rank]
            summed_likelihood = df_l["log_likelihood"].sum()
            data = {rank_by:rank,
                    "Root":root,
                    "summed_likelihood":summed_likelihood}
            df_totals = df_totals.append(data, True)

            total_rec = count*number_of_roots
            pct_progress =  round((progress/total_rec)*100,2)

            if (rank % LOG_EVERY_N) == 0:
                print("Current root",root,", Overall progress = %",pct_progress)
            progress += 1

    timestr = time.strftime("%d%m%y")
    file_name = rank_by + "_" + timestr + "_summed.csv"
    df_totals.to_csv(file_name)

    print("\nSummed_likelihoods have been written to:",file_name)


    print("\n*******************************************","\n")
    return df_totals

def summed_likelihood_high_rep(rank_by,species_representation,count):

    """ Function to sum likelihoods based on rank for each root """

    print("\n*******************************************","\n")
    print("\nSummed likelihoods and ranks now being calculated for high" + \
          " species representation uml_rec files")

    LOG_EVERY_N = 100
    number_of_roots = 0
    rtt  = open("roots_to_test.txt", "r")
    for root in rtt:
        number_of_roots += 1
    rtt.close()
    print("Number of roots to sum likelihood for:", number_of_roots)

    timestr = time.strftime("%d%m%y")
    file_name = rank_by + "_" + timestr + "_high_rep_ranked.csv"
    df = pd.read_csv(file_name)
    d = {rank_by: [], 'Root': [], "summed_likelihood": [] }
    df_totals = pd.DataFrame(data=d)

    rtt  = open("roots_to_test.txt", "r")
    progress = 0
    for root in rtt:
        for rank in range(0,len(species_representation),1):
            root = root.strip()
            df_r = df[df["Root"] == root]
            df_l = df_r[df_r[rank_by]>= rank]
            summed_likelihood = df_l["log_likelihood"].sum()
            data = {rank_by:rank,
                    "Root":root,
                    "summed_likelihood":summed_likelihood}
            df_totals = df_totals.append(data, True)

            total_rec = count*number_of_roots
            pct_progress =  round((progress/total_rec)*100,2)

            if (rank % LOG_EVERY_N) == 0:
                print("Current root",root,", Overall progress = %",pct_progress)
            progress += 1

    timestr = time.strftime("%d%m%y")
    file_name = rank_by + "_" + timestr + "_high_rep_summed.csv"
    df_totals.to_csv(file_name)

    print("\nSummed_likelihoods for high rep files have been written to:",
          file_name)
    print("\n*******************************************","\n")
    return df_totals

def plot_results():

    """ Function to plot results """

    print("\n##### Results are being plotted #####")
    timestr = time.strftime("%d%m%y")
    file_name = rank_by + "_" + timestr + "_summed.csv"
    df = pd.read_csv(file_name)
    df = df.sort_values(by=rank_by, ascending=False)
    df["ML_root"] = np.where(df['Root']==ML_root, 'yes', 'no')
    order = list(df[rank_by])

    sns.set_style("whitegrid", {'axes.grid' : False})
    sns.set_palette("flare")

    fig = plt.figure()
    fig.set_size_inches(15, 10)
    ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
    ax = sns.lineplot(data=df, x=rank_by, y="summed_likelihood",
                      palette="flare",
                      hue = "Root",
                      style = "ML_root",
                      style_order = ["no","yes"])
    ax.set_xlabel(rank_by + "_rank", size = 12)
    ax.set_ylabel('summed_likelihood', size = 12)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    fig_name = rank_by + "_" + timestr + "_summed.png"
    ax.figure.savefig(fig_name)

    print("\n##### High species representation results are being plotted #####")
    timestr = time.strftime("%d%m%y")
    file_name = rank_by + "_" + timestr + "_high_rep_summed.csv"
    df = pd.read_csv(file_name)
    df = df.sort_values(by=rank_by, ascending=False)
    df["ML_root"] = np.where(df['Root']==ML_root, 'yes', 'no')
    order = list(df[rank_by])

    sns.set_style("whitegrid", {'axes.grid' : False})
    sns.set_palette("flare")

    fig = plt.figure()
    fig.set_size_inches(15, 10)
    ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
    ax = sns.lineplot(data=df, x=rank_by, y="summed_likelihood",
                      palette="flare",
                      hue = "Root",
                      style = "ML_root",
                      style_order = ["no","yes"])
    ax.set_xlabel(rank_by + "_rank", size = 12)
    ax.set_ylabel('summed_likelihood', size = 12)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    fig_name = rank_by + "_" + timestr + "_high_species_summed.png"
    ax.figure.savefig(fig_name)

def clean_up(rank_by):
    print("\n*******************************************","\n")

    file_name = rank_by + "_ratio_results"
    os.mkdir(file_name)
    timestr = time.strftime("%d%m%y")
    command = "mv *" + rank_by + "* " + file_name
    os.system(command)


def output():
    print("\n*******************************************","\n")
    print("Software has finished running \n")
    print("All figures are saved as *.png " + \
          "and all data has been saved as *.csv files to the ratio_results dir")

def main():
    usage()
    count = get_count(ML_root)
    species_representation = get_species_rep(ML_root)
    df = run_DTL(count,ML_root,rank_by,species_representation)
    write_results(df,rank_by)
    make_high_rep_df(df,rank_by)
    df_totals = summed_likelihood(rank_by,count)
    df_total_high_rep = summed_likelihood_high_rep(rank_by,
                                                   species_representation,count)
    plot_results()
    clean_up(rank_by)
    output()

if __name__ == "__main__":
    ML_root = sys.argv[1]
    rank_by = sys.argv[2]
    main()
