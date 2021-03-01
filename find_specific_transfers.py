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

def get_transfers(uTs_file):

    """ Function to find frequency of each transfer """

    t_from_list = []
    t_to_list = []
    t_freq_list = []


    uTs_file_open = open(uTs_file, "r")
    for line in uTs_file_open:
        if line.startswith("#"):
            continue
        else:
            line = line.strip()
            line_list = line.split("\t")

            t_from = line_list[0]
            t_to = line_list[1]
            t_freq = float(line_list[2])

            t_from_list.append(t_from)
            t_to_list.append(t_to)
            t_freq_list.append(t_freq)

    t_file_list = [uTs_file]*len(t_from_list)

    data = {"Gene_family":t_file_list, "Transfer_from":t_from_list,
            "Transfer_to":t_to_list, "Transfer_frequency":t_freq_list}
    df = pd.DataFrame.from_dict(data)

    return df

def find_likely_transfers(dataframe):

    """ Function to find the most likely transfers """

    df_likely = dataframe[dataframe["Transfer_frequency"] >= 0.80]

    return(df_likely)

def run_transfer_finder():

    """ Function to execute DTL calculation functions """

    print("\n*****************************************************************")
    print("\nCalculating the transfer frequencies for all transfer")

    data = {"Gene_family":[], "Transfer_from":[],
            "Transfer_to":[], "Transfer_frequency":[]}

    df_all = pd.DataFrame(data)
    for file_name in glob.glob("*.uTs"):
        df = get_transfers(file_name)
        df_all = df_all.append(df, True)
    df_all.to_csv("All_transfer_frequencies.csv")

    print("\n*****************************************************************")
    print("\nCalculating the most likely transfers")

    df_likely = find_likely_transfers(df_all)
    df_likely.to_csv("likely_transfers.csv")

def output():

    """ Function to explain output"""

    print("\n*****************************************************************")
    print("\nAll transfer data is saved in 'All_transfer_frequencies.csv'")
    print("\nThe most likely transfers are saved in 'likely_transfers.csv'")

    os.mkdir("Transfer_results")
    os.system("mv *.csv Transfer_results")

    print("\nBoth results are saved in the 'Transfer_results' directory")
    print("\nScript finished running")
    print("\n*****************************************************************")

def main():
    run_transfer_finder()
    output()

if __name__ == "__main__":
    main()
