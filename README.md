# ALE_methods #

* Python util scripts to facilitate analysis of duplication, transfer and loss with Amalgamated Likelihood Estimation (ALE). 

## Making ALE objects ##
1. Genomes are downloaded from repositories in FASTA format.
   1. Genome are reformatted to be compatible with ALE 
   2. Same genomes in this fashion: “ecoli.fa” for genus only or “Arabidopisis.thaliana.fa” for genus and species.
   3. Use _rf_genome_for_mcl.py_ to reformat genomes for ALE.
        - Code requires one command line argument: a species list in the form of a text file e.g. “species_list_demo.txt"

1. Construct pairwise alignments
   1. Perform all vs all Diamond BLAST https://github.com/bbuchfink/diamond
   2. Merge all the genomes into a single fasta file with the following Bash command: _cat *.fa > ALE_demo_dataset.fasta_
   3. Make a diamond blast database. Bash command: _diamond makedb --in ALE_demo_dataset.fasta -d nr_
   4. Blast ALE_demo_dataset.fasta against diamond blast database with the bash command (Adjust parameters to suit data and machine): _diamond blastp -d nr -q ALE_demo_dataset.fasta -o AvA_demo.txt -f 6 -p 40 –evalue 0.0001 –log_
        - Output is in AvA_demo.txt

1. Cluster gene families using MCL https://micans.org/mcl/.
   1. Extract the results required from AvA_demo.txt using the bash command: _cut -f 1,2,11 AvA_demo > seq.abc_
   2. Use MCL to first load the blast results, with the bash command: _mcxload -abc seq.abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(2002)' -o seq.mci -write-tab seq.tab_


