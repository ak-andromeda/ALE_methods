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


def get_DTL(rec_file):

    """ Function to extract the loss/speciation ratio from rec_file """

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

    return t

def find_transfers():

    """ Function to find uml_rec with transfers """

    print("\n*****************************************************************")
    print("\nTransfer rates are being obtained")

    transfer_dict = {}

    for file in glob.glob("*.uml_rec"):
        rec_file = open(file, "r")
        transfer_rate = get_DTL(rec_file)
        transfer_dict[file] = transfer_rate

    return transfer_dict

def make_dataframe(transfer_dict):

    """ Function to construct dataframe from transfer dict """

    print("\n*****************************************************************")
    print("\nTransfer rates are analysed")

    rec_files = []
    transfer_rates = []

    for key, val in transfer_dict.items():
        rec_files.append(key)
        transfer_rates.append(val)

    data = {"rec_file":rec_files,"transfer_rates":transfer_rates}
    df = pd.DataFrame.from_dict(data)
    df.to_csv("transfer_rates.csv")


def main():
    transfer_dict = find_transfers()
    make_dataframe(transfer_dict)

if __name__ == "__main__":
    main()
