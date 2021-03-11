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
    file_name = "Species_repsentation_" + timestr + ".txt"
    file = open(file_name, "w")
    for key, val in species_rep_dictionary.items():
        line = key + "," + str(val) + "\n"
        file.write(line)
    file.close()

    print("\nSpecies representation calculated and saved in the text file:",
           file_name)

    print("\n*******************************************","\n")

    return(species_rep_dictionary)

def get_log_likelihood(rec_file):

    """ Function to extract the loss/speciation ratio from rec_file """

    for l in rec_file:
        if l.startswith(">logl:"):
            results = l.strip()
            results = float(results[7:])

    return results

def sort_dict(dict):

    """ Function to sort dictionary for lowest to highest """

    sorted_dict = {k:v for k,v in sorted(dict.items(),key=lambda item: item[1])}
    return(sorted_dict)

def get_rec_file_rank(sorted_values):

    """ Function to assign rank to orthogroup based on metric"""

    number_rec_files = 0
    rec_file_rank = {}
    for file_name, metric in sorted_values.items():
        rec_file_rank[file_name] = number_rec_files
        number_rec_files +=1

    return rec_file_rank


def get_rank_for_all(count,ML_root,ranked_rec_file,species_representation):

    """ function to execute DTL calculation functions """

    print("DTL for each gene family are being calculated")

    rtt  = open("roots_to_test.txt", "r")

    d = {"Root": [], "Gene_family": [], "log_likelihood": [],
         "species_rep_rank": [] }

    df = pd.DataFrame(data=d)
    colums = list(df)


    root_num = 0
    for root in rtt:
        root = root.strip()
        root_num += 1
        progress = 0
        print("\nRoot currently under analysis:", root)
        new_path = os.getcwd() + "/" + root
        with cd(new_path):
            current_dir = os.getcwd()
            print(current_dir)
            for file_name in glob.glob("*.uml_rec"):
                rec_file = open(file_name, "r")
                log_likelihood = get_log_likelihood(rec_file)
                species_present = species_representation[file_name]
                rank = ranked_rec_file[file_name]
                data = {"Root":root,
                        "Gene_family":file_name,
                        "species_representation": species_present,
                        "log_likelihood":log_likelihood,
                        "species_rep_rank":rank}
                df = df.append(data, True)
                rec_file.close()
                progress +=1
                print("root_num",root_num,"progress",progress,"/",count)
    return(df)

def write_results(df):

    """ Function to write pandas to csv """

    timestr = time.strftime("%d%m%y")
    file_name = "Species_rep_" + timestr + "_ranked.csv"
    df.to_csv(file_name)
    print("\nResults have been written to",  file_name)
    print("\nSlice of dataframe shown below:")
    print(df.head())

def summed_likelihood(count):

    """ Function to sum likelihoods based on rank for each root """

    print("\n*******************************************","\n")
    print("\nSummed likelihoods and ranks now being calculated" + \
          " for uml_rec files\n")

    timestr = time.strftime("%d%m%y")
    file_name = "Species_rep_" + timestr + "_ranked.csv"
    df = pd.read_csv(file_name)
    d = {'Root': [], "species_rep_rank": [],
         "summed_likelihood": [], "ML_diff":[] }
    df_totals = pd.DataFrame(data=d)

    rtt  = open("roots_to_test.txt", "r")
    root_num = 0
    for root in rtt:
        root_num += 1
        progress = 0
        for rank in range(0,int(count),1):
            progress += 1
            root = root.strip()
            df_r = df[df["Root"] == root]
            df_l = df_r[df_r["species_rep_rank"] >= rank]
            summed_likelihood = df_l["log_likelihood"].sum()

            df_ml = df[df["Root"] == ML_root]
            df_mlr = df_ml[df_ml["species_rep_rank"] >= rank]
            summed_likelihood_ml = df_mlr["log_likelihood"].sum()
            ML_diff = summed_likelihood - summed_likelihood_ml
            data = {"Root":root,
                    "species_rep_rank":rank,
                    "summed_likelihood":summed_likelihood,
                    "ML_diff":ML_diff}

            df_totals = df_totals.append(data, True)
            print("root_num",root_num,"progress",progress,"/",count)


    timestr = time.strftime("%d%m%y")
    file_name = "Species_rep_" + timestr + "_summed.csv"
    df_totals.to_csv(file_name)

    print("\nSummed_likelihoods have been written to:",file_name)


    print("\n*******************************************","\n")
    return df_totals

def plot_results():

    """ Function to plot results """

    print("\n##### Results are being plotted #####")
    timestr = time.strftime("%d%m%y")
    file_name = "Species_rep_" + timestr + "_summed.csv"
    df = pd.read_csv(file_name)
    df = df.sort_values(by="species_rep_rank", ascending=True)
    df["ML_root"] = np.where(df['Root']==ML_root, 'yes', 'no')

    sns.set_style("whitegrid", {'axes.grid' : False})
    sns.set_palette("flare")

    fig = plt.figure()
    fig.set_size_inches(15, 10)
    ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
    ax = sns.lineplot(data=df, x="species_rep_rank", y="summed_likelihood",
                      palette="husl",
                      hue = "Root",
                      style = "ML_root",
                      style_order = ["no","yes"])
    ax.set_xlabel("Species_rep_rank", size = 12)
    ax.set_ylabel('summed_likelihood', size = 12)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    fig_name ="Species_rep_" + timestr + "_summed.png"
    ax.figure.savefig(fig_name)

    fig = plt.figure()
    fig.set_size_inches(15, 10)
    ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
    ax = sns.lineplot(data=df, x="species_rep_rank", y="ML_diff",
                      palette="husl",
                      hue = "Root",
                      style = "ML_root",
                      style_order = ["no","yes"])
    ax.set_xlabel("Species_rep_rank", size = 12)
    ax.set_ylabel('Summed_likelihood_diff_from_ML', size = 12)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    fig_name ="Species_rep_ml_diff" + timestr + "_summed.png"
    ax.figure.savefig(fig_name)

def clean_up():
    print("\n*******************************************","\n")

    file_name = "species_rep_ratio_results"
    os.mkdir(file_name)
    timestr = time.strftime("%d%m%y")
    command = "mv *Species_rep* " + file_name
    os.system(command)

def main():
    count = get_count(ML_root)
    #species_representation = get_species_rep(ML_root)
    #sorted_species_representation = sort_dict(species_representation)
    #ranked_rec_file = get_rec_file_rank(sorted_species_representation)
    #df = get_rank_for_all(count,ML_root,ranked_rec_file,species_representation)
    #df = write_results(df)
    #df_totals = summed_likelihood(count)
    #plot_results()
    clean_up()

if __name__ == "__main__":
    ML_root = sys.argv[1]
    main()
