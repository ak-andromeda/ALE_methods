import sys
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import time
import glob
import Bio
from Bio import SeqIO

def calculate_size_rep(ALE_species_list):

	""" Function to calculate size and species represention
	 	requires a species list and grouping text file e.g.
		< Ecoli,bacteria > """

	species_count_total = {}
	species_list = []
	species_group_dict = {}
	file = open(ALE_species_list, "r")
	for l in file:
		l = l.strip()
		temp_list = l.split(",")
		species = temp_list[0]
		group = temp_list[1]
		species_list.append(species)
		species_count_total[species] = 0
		species_group_dict[species] = group

	orthogroups = []
	size = []
	size_species_rep = []
	for orthogroup in glob.glob("*.aln"):
		orthogroups.append(orthogroup)
		records = SeqIO.parse(orthogroup, "fasta")
		number_records = 0

		species_dict = {}
		for record in records:
			number_records += 1
			record_desc = record.description
			for species in species_list:
				if species in record_desc:
					species_dict[species] = 1
					species_count_total[species] +=1
				else:
					continue

		species_rep = (len(species_dict)/len(species_list))*100
		size_species_rep.append(species_rep)
		size.append(number_records)

	species = []
	total_records = []
	groups = []
	for key,val in species_count_total.items():
		species.append(key)
		total_records.append(val)
		group = species_group_dict[key]
		groups.append(group)

	data_totals = {'Species': species, "Total_records": total_records,
				   'Group':groups}

	df_totals = pd.DataFrame.from_dict(data_totals)
	print(df_totals.head())
	print(df_totals.tail())
	df_totals.to_csv("Species_totals.csv")
	print("\n")

	data = {'Orthogroup':orthogroups,'Size':size,'Species_rep':size_species_rep}
	df = pd.DataFrame.from_dict(data)
	df.to_csv("Orthogroup_size_rep.csv")
	print(df.head())

def plot_dist(dataset_name, paper_pal):

	""" Function to plot Distribution """

	timestr = time.strftime("%d%m%y")

	df = pd.read_csv("Orthogroup_size_rep.csv")
	print(df.head())

	sns.set_style("whitegrid", {'axes.grid' : False})
	g = sns.displot(data=df.Size, color = paper_pal[0])
	g.set_axis_labels("Gene family size", "Count", labelpad=12)
	g.fig.set_size_inches(6, 6)
	fig_name = dataset_name + "_Hist_" + timestr + ".png"
	g.savefig(fig_name)

def plot_kde(dataset_name, paper_pal):

	""" Function to plot Distribution """

	timestr = time.strftime("%d%m%y")

	df = pd.read_csv("Orthogroup_size_rep.csv")
	print(df.head())

	sns.set_style("whitegrid", {'axes.grid' : False})
	g = sns.displot(data=df, x="Size", y="Species_rep", kind = "kde",
					fill=True, color = paper_pal[0])
	g.set_axis_labels("Gene family size", "Species represention", labelpad=12)
	g.fig.set_size_inches(6, 6)
	fig_name = dataset_name + "_KDE_" + timestr + ".png"
	g.savefig(fig_name)

def plot_lm(dataset_name, paper_pal):

	""" Function to plot Distribution """

	timestr = time.strftime("%d%m%y")

	df = pd.read_csv("Orthogroup_size_rep.csv")
	print(df.head())

	sns.set_style("whitegrid", {'axes.grid' : False})
	sns.set_palette(paper_pal)
	g = sns.lmplot(data=df, x="Size", y="Species_rep", palette = paper_pal,
				   line_kws={"color": paper_pal[1]})
	g.set_axis_labels("Gene family size", "Species represention", labelpad=12)
	g.fig.set_size_inches(6, 6)
	g.set(ylim=(0, 110))
	fig_name = dataset_name + "_LM_" + timestr + ".png"
	g.savefig(fig_name)

def plot_totals(dataset_name,paper_pal):

	""" Function to plot Distribution """

	timestr = time.strftime("%d%m%y")

	df = pd.read_csv("Species_totals.csv")
	print(df.head())

	sns.set_style("whitegrid", {'axes.grid' : False})
	g = sns.displot(df, x="Total_records", hue="Group", multiple="stack",
						kind = "kde", legend = True,
						palette = paper_pal)

	g.set_axis_labels("Genome size", "Density", labelpad=12)
	g.fig.set_size_inches(12, 6)
	fig_name = dataset_name + "_species_totals_" + timestr + ".png"
	g.savefig(fig_name)

def plot_subplot(dataset_name,paper_pal):

	""" Function to plot subplot """

	timestr = time.strftime("%d%m%y")

	df = pd.read_csv("Orthogroup_size_rep.csv")

	sns.set_style("whitegrid", {'axes.grid' : False})

	# Set up the matplotlib figure
	f, axes = plt.subplots(2, 2, figsize=(12,12), sharex=False)

	# Plot a simple histogram with binsize determined automatically
	sns.histplot(data = df.Size, bins=20, color=paper_pal[1], ax=axes[0,0])

	# Plot a simple histogram with binsize determined automatically
	sns.histplot(data = df.Species_rep, bins=20, color=paper_pal[1], ax=axes[0,1])

	# Plot a kernel density estimate
	g = sns.kdeplot(data=df, x="Size", y="Species_rep", kind = "kde",
			        fill=True, color=paper_pal[0], ax = axes[1,0])
	g.set(xlim=(0,100))


	# Plot a lm plot
	g = sns.regplot(data=df, x="Size", y="Species_rep", color = paper_pal[0],
				line_kws={"color":paper_pal[1]}, ax=axes[1,1])
	g.set(ylim=(0, 120))

	plt.savefig('subplot_distributions.png')

def main():
	#paper_pal = ["#feaea5","#fe8976","#fc986e","#fcb29c",
				 #"#a8fbf8","#b1bffb","#c4c4c4"]
	paper_pal = ["#fc986e","#c4c4c4"]
	calculate_size_rep(species_list)
	plot_dist(dataset_name,paper_pal)
	plot_kde(dataset_name,paper_pal)
	plot_lm(dataset_name,paper_pal)
	plot_totals(dataset_name, paper_pal)
	plot_subplot(dataset_name,paper_pal)

if __name__ == "__main__":
	species_list = sys.argv[1]
	dataset_name = sys.argv[2]
	main()
