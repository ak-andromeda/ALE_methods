import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def plot_subplot(df,paper_pal):

    """ Function to plot subplot """

    sns.set_style("whitegrid", {"axes.grid":False})

    # Set up the matplotlib figure
    f, axes = plt.subplots(3, 1, figsize=(12,12), sharex=False)

    sns.histplot(data = df.Fraction_missing,
                 color=paper_pal[1], ax=axes[0])
    sns.kdeplot(data=df, x="Fraction_missing", hue="Group",
                fill=True, common_norm=False, palette="flare",
                alpha=.65, linewidth=0, ax=axes[1])
    sns.kdeplot(data=df, x="Fraction_missing", hue="Clade",
                fill=True, common_norm=False, palette="flare",
                alpha=.75, linewidth=0,
                multiple="stack", ax=axes[2])

    plt.savefig('Fraction_missing_subplot_distributions.png')

def plot_bar(df,paper_pal):

    """ Plot a bar plot of fraction missing """

    sns.set_style("whitegrid", {"axes.grid":False})
    sns.set(rc={'figure.figsize':(12,14)})
    sns.barplot(x="Species",y="Fraction_missing",data = df,
                palette = "flare",color=paper_pal[1])
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig('Fraction_missing_species_bar.png')


def make_df(species_list,group_list,clade_list):

     """Function to make df from list"""

     group_dict = {}
     file = open(group_list,"r")
     for line in file:
         line = line.strip().split(",")
         species = line[0]
         group = line[1]
         group_dict[species] = group

     clade_dict = {}
     file = open(clade_list,"r")
     for line in file:
         line = line.strip().split(",")
         species = line[0]
         clade = line[1]
         clade_dict[species] = clade

     file = open(species_list, "r")
     species_list = []
     fraction_missing_list = []
     groups = []
     clades = []
     for line in file:
         line = line.strip().split(":")
         species = line[0]
         fm = line[1]
         group = group_dict[species]
         clade = clade_dict[species]
         species_list.append(species)
         fraction_missing_list.append(float(fm))
         clades.append(clade)
         groups.append(group)
     file.close()

     data = {"Species":species_list,
             "Fraction_missing":fraction_missing_list,
             "Group":groups,
             "Clade":clades}
     df = pd.DataFrame.from_dict(data)
     df = df.sort_values(by=['Fraction_missing'])
     return df


def main():
    paper_pal = ["#fc986e","#c4c4c4", "#feaea5","#fe8976","#fcb29c",
                 "#a8fbf8","#b1bffb"]
    df = make_df(species_list,group_list,clade_list)
    plot_bar(df, paper_pal)
    plot_subplot(df, paper_pal)




if __name__ == "__main__":
    species_list = sys.argv[1]
    group_list = sys.argv[2]
    clade_list = sys.argv[3]
    main()
