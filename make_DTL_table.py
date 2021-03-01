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

def get_DTL(rec_file):

    """ Function to extract the loss/speciation ratio from rec_file """

    DTLS_return = []
    log_l = []

    for l in rec_file:
        if l.startswith("Total"):
            results = l.strip()
            result_list = results.split(" ")
            DTLS = result_list[1]
            DTLS_list = DTLS.split("\t")[1:]

            DTLS_return.append(float(DTLS_list[0]))
            DTLS_return.append(float(DTLS_list[1]))
            DTLS_return.append(float(DTLS_list[2]))
            DTLS_return.append(float(DTLS_list[3]))

        elif l.startswith(">logl:"):
            log_l.append(float(l[7:]))

    return DTLS_return, log_l[0]

def run_DTL():

    """ function to execute DTL calculation functions """

    print("\n*****************************************************************")
    print("\nRunning DTL")

    d = {'Root':[], 'Gene_family':[], "log_likelihood":[],
         'duplication':[], 'transfer':[], 'loss':[], 'speciation':[]}

    df = pd.DataFrame(data=d)
    colums = list(df)
    rtt  = open("roots_to_test.txt", "r")
    for root in rtt:
        root = root.strip()
        print("\nRoot currently under analysis:", root)
        new_path = os.getcwd() + "/" + root
        with cd(new_path):
            current_dir = os.getcwd()
            print(current_dir)
            count = 0
            for file_name in glob.glob("*.uml_rec"):
                count += 1
                print(count,"/",30625)
                rec_file = open(file_name, "r")
                DTLS, logl = get_DTL(rec_file)

                data = {"Root":root,
                        "Gene_family":file_name,
                        "log_likelihood":logl,
                        "duplication":DTLS[0],
                        "transfer":DTLS[1],
                        "loss":DTLS[2],
                        "speciation":DTLS[3]}

                df = df.append(data, True)
                rec_file.close()
            print("Root finished")

            print(df.head())
            print(df.tail())

    df.to_csv("DTL_and_likelihoods_table.csv")
    return(df)

def main():
    df = run_DTL()
    print(df.head())
    print(df.tail())

if __name__ == "__main__":
    main()
