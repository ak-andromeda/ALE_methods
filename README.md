# ALE_methods #

* Python util scripts to facilitate analysis of duplication, transfer and loss with Amalgamated Likelihood Estimation (ALE). 

## Making ALE objects ##
1. Genomes are downloaded from repositories in FASTA format.
   1. Genome are reformatted to be compatible with ALE 
   2. Same genomes in this fashion: “ecoli.fa” for genus only or “Arabidopisis.thaliana.fa” for genus and species.
   3. Use _rf_genome_for_mcl.py_ to reformat genomes for ALE.
      - Code requires one command line argument: a species list in the form of a text file e.g. “species_list_demo.txt"

1. Construct pairwise alignments
   1. Perform all vs all Diamond BLAST (https://github.com/bbuchfink/diamond)
   2. Merge all the genomes into a single fasta file with the following bash command: _cat *.fa > ALE_demo_dataset.fasta_
   3. Make a diamond blast database. Bash command: _diamond makedb --in ALE_demo_dataset.fasta -d nr_
   4. Blast ALE_demo_dataset.fasta against diamond blast database with the bash command (Adjust parameters to suit data and machine): _diamond blastp -d nr -q ALE_demo_dataset.fasta -o AvA_demo.txt -f 6 -p 40 –evalue 0.0001 –log_
      - Output is in AvA_demo.txt

1. Cluster gene families using MCL (https://micans.org/mcl/).
   1. Extract the results required from AvA_demo.txt using the bash command: _cut -f 1,2,11 AvA_demo > seq.abc_
   2. Use MCL to first load the blast results, with the bash command: _mcxload -abc seq.abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(2002)' -o seq.mci -write-tab seq.tab_
   3. Use MCL to cluster the blast results with the folliwng bash command: _mcl seq.mci -I 1.4 -use-tab seq.tab._
      - Ensure MCL is added to path.
      - Step 3. can be repeated using different inlfation values (-I parameter) to adjust gene family cluster sizes.
      - Output is in a single text file: 'out.seq.mci.I14'
   4. Convert MCL output file to individual FASTA files.
      - _Use make_gene_families_mcl.py_ to convert the MCL output into individual FASTA files.
      - The script requires two command line arguments:  1: the output from the MCL ('out.seq.mci.I14'), 2: The dataset used to in the all vs all blast ('ALE_demo_dataset.fasta').

1. Align and trim gene families
   1. Use orthoprocess.py to remove families with less than four sequences present. 
      - Four is the minimum number of sequences in which MAFFT can align - removing these files early helps audit the process.
   2. Use MAFFT (https://mafft.cbrc.jp/alignment/software/) to align the gene families. 
      - This can be done with GNU parallel (https://www.gnu.org/software/parallel/ to speed up the process).
      - The following code can be used to align in parallel: _ls *.fasta | parallel mafft –auto {} “>” {}.aln._
   3. Use BMGE (ftp://ftp.pasteur.fr/pub/GenSoft/projects/BMGE/) to trim poorly aligned sites. 

1. Infer bootstrap distribution of trees for each gene family. This can be done with IQ-Tree (http://www.iqtree.org/).
   1. For protein sequences, the LG in combination with the C10:C60* provides a parameter rich model to infer the bootstrap distribution of trees. Please see IQ-Tree documentation for more help.
      - bash command to run IQ-Tree: i_qtree -s <gene_family_ID.aln> -m MFP -mset LG+C20,LG+C30,LG+C40 -madd -bb 1000 -wbtl 
      - Additional paramaters can be added such as 'LG+C20+g+f' to model for site rate heterogeneity.  


