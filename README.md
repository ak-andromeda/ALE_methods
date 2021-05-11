# ALE_methods #

* Python util scripts to facilitate analysis of duplication, transfer and loss with Amalgamated Likelihood Estimation (ALE). 

## Make ALE objects ##
1. Genomes are downloaded from repositories in FASTA format.
   1. Genome are reformatted to be compatible with ALE 
   2. Name genomes in this fashion: “Sulfolobus.fa” for genus only or “Arabidopisis.thaliana.fa” for genus and species.
   3. Use _rf_genome_for_mcl.py_ to reformat genomes for ALE.
      - Code requires one command line argument: a species list in the form of a text file e.g. “species_list_demo.txt"
      - Genome names must match species list. e.g. if species list includes 'Ecoli', genome in the directory must be named 'Ecoli.fa'. 

1. Construct pairwise alignments
   1. Perform all vs all Diamond BLAST (https://github.com/bbuchfink/diamond)
   2. Merge all the genomes into a single fasta file with the following bash command: _cat *_rf.fa > ALE_demo_dataset.fasta_
   3. Make a diamond blast database. Bash command: _diamond makedb --in ALE_demo_dataset.fasta -d nr_
   4. Blast ALE_demo_dataset.fasta against diamond blast database with the bash command (Adjust parameters to suit data and machine): _diamond blastp -d nr -q ALE_demo_dataset.fasta -o AvA_demo.txt -f 6 -p 40 –evalue 0.0001 –log_
      - Output is in AvA_demo.txt

1. Cluster gene families using MCL (https://micans.org/mcl/).
   1. Extract the results required from AvA_demo.txt using the bash command: _cut -f 1,2,11 AvA_demo.txt > seq.abc_
   2. Use MCL to first load the blast results, with the bash command: _mcxload -abc seq.abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(2002)' -o seq.mci -write-tab seq.tab_
   3. Use MCL to cluster the blast results with the following bash command: _mcl seq.mci -I 1.4 -use-tab seq.tab._
      - Ensure MCL is added to path.
      - Step 3. can be repeated using different infltion values (-I parameter) to adjust gene family cluster sizes.
      - Output is in a single text file: 'out.seq.mci.I14'
   4. Convert MCL output file to individual FASTA files.
      - _Use make_gene_families_mcl.py_ to convert the MCL output into individual FASTA files.
      - The script requires two command line arguments:  1: the output from the MCL ('out.seq.mci.I14'), 2: The dataset used to in the all vs all blast ('ALE_demo_dataset.fasta').

1. Align and trim gene families
   1. Use _find_viable_gene_families.py_ to remove families with less than four sequences present. 
      - Four is the minimum number of sequences in which MAFFT can align - removing these files early helps audit the process.
   2. Use MAFFT (https://mafft.cbrc.jp/alignment/software/) to align the gene families. 
      - This can be done with GNU parallel (https://www.gnu.org/software/parallel/ to speed up the process).
      - First change into the directory with the viable sequences _cd four_plus_seq/_
      - The following code can be used to align in parallel: _ls *.fasta | parallel mafft –auto {} “>” {}.aln._
   3. Use BMGE (ftp://ftp.pasteur.fr/pub/GenSoft/projects/BMGE/) to trim poorly aligned sites. 

1. Infer bootstrap distribution of trees for each gene family. This can be done with IQ-Tree (http://www.iqtree.org/).
   1. For protein sequences, the LG in combination with the C10:C60* provides a parameter rich model to infer the bootstrap distribution of trees. Please see IQ-Tree documentation for more help.
      - bash command to run IQ-Tree: _iqtree -s <gene_family_ID.aln> -m MFP -mset LG+C20,LG+C30,LG+C40 -madd -bb 1000 -wbtl_
      - Additional paramaters can be added such as 'LG+C20+g+f' to model for site rate heterogeneity.  
1. Convert bootstrap tree distributions to ALE object using ALEobserve.
   1. To create ALE object use: _ALEobserve < ufboot file >_
   2. To speed up the process use GNU parallel: _ls *.ufboot | parallel ALEobserve {}.ale_

## Infer unrooted species tree ##

1. Use OrthoFinder (https://github.com/davidemms/OrthoFinder) to infer single copy orthologs.
   1. See OrthoFinder documentation for instructions. 
   2. If no single copy orthologs were identified (this is often the case for multicellular organisms). Use _prem3.py_ to generate single copy gene families from the orthologous sequence files produced by Orthofinder. 
      - _prem3.py_ takes two command line arguments: 1) A species list ('species_list_demo.txt'), 2) A percentage cutoff for species representation remaining of each gene family. 
      - Adjust the percentage cutoff to determine the number of single copy families returned. Aim to increase the value as close to '100' as possible, whilst retaining a sufficient amount of data to construct a sequence alignment. 

1. Construct a super matrix alignment, or employ a super tree approach. 
   1. Align single copy gene families with MAFFT (As shown in the Make ALE Objects section). 
   2. Trim gene families with BMGE to remove poorly aligned sites. 
   3. Use _super_matrix_2.py_ to construct a super matrix alignment.
      - _super_matrix_2.py_ takes one command line argument, a species list e.g. 'species_list_demo.txt'
      - The species list must match with the Fasta headers of the sequences in the gene families. 
      - This will be the case of all the genomes are formatted as above. 
   4. Infer an unrooted species tree:
      - Use IQ-Tree as described above to infer a species tree. 
      - Alternatively, instead of constructing a super matrix, infer individual gene trees for the single copy gene families and employ a super tree approach such as ASTRAL (https://github.com/smirarab/ASTRAL) 
   5. Root the tree in different candidate position. 
      - To visualise unrooted trees to help select candidate root positions use ITOL (https://itol.embl.de/) or FigTree (http://tree.bio.ed.ac.uk/software/figtree/)

## Reconcile ALE objects with candidate rooted species trees ##

1. Reconcile ale objects with each candidate rooted species tree.
   1. Use ALEml_undated to reconcile each ALE object with tge candidate rooted species trees
      - Use bash command: ALEml_undated <rooted_species_tree> <ale_object>
      - You perform this command for every ale_object, under the different candidate root positions. 
      - e.g bash command to run ALE: _ALEml_undated root_on_ecoli.newick 1.ale_
      - You can also factor in the percentage of the genome that is missing using BUSCO (https://busco.ezlab.org/)
      - BUSCO provides a completion value e.g 80% - thus, the fraction missing is 0.20. See 'fraction_missing.txt' for e.g.
      - To add the fraction missing missing parameter just add: _fraction_missing=fraction_missing.txt_ (where 'fraction_missing.txt' holds the respective fraction missing for each sequence e.g. 'Ecoli:0.20'

1. Perform an Approximately unbiased (AU) test (Shimodiara, 2002) to identify likely root position.
   1. There should be three output files from ALEml_undated: ‘.uml_rec’, ‘.uTs’ and ‘.tree’.
      - The ‘.uml_rec’ files generated for each candidate root should be in separate directories. 
      - The naming of each ‘.uml_rec’ file should be consistent for each root. I.e. if you have tested 4 roots, in directories named root_1, root_2, root_3 and root_4, the reconciliation output for the ALE object: 1.ale, should be named 1.ale.uml_rec in each of the four directories. 1.ale.uml_rec will have different values for DTL under each candidate root.
   2. Use the write_consel_file_p3.py script to construct a table of likelihoods. 
      - The script requires the list of directories holding the ‘.uml_rec’ for each candidate root position. Save the output of the script into an output file.
      - For e.g.: _write_consel_file_p3.py root_1/ root_2/ root_3/ root_4/ > likelihoods_table_
      - The order of the roots is not preserved in the likelihood table. Open the likelihood table in a text editor and note the order of the roots in the likelihood table. 
      - A python2 version is also available _write_conse_file.py_
   3. Use the following bash command to copy the table into a .mt file: _cp likelihoods_table likelihoods_table.mt _
   4. Use Consel (http://stat.sys.i.kyoto-u.ac.jp/prog/consel/) to perform an AU test to assess the likelihoods of each root. 
      - Consel requires three commands to be executed in succession 
      - _consel/bin/makermt likelihoods_table_
      - _consel/bin/consel likelihoods_table_
      - _consel/bin/catpv likelihoods_table > au_test_out_
   5. Assess the likelihood of each root in the 'au_test_out' file 
      - The ‘au_test_out.txt’ file includes a table ranking the roots by likelihood, and a number of statistical tests. The AU value column provides a ‘p-value’. If the AU p-value is <0.05 the candidate roots position has been significantly rejected. Roots that have a p-value greater than 0.05 can not be significantly rejected and are therefore likely root positions for the phylogenetic tree. See 'au_test_out.txt' as an example.

## Post-hoc alaysis ##

1. Assessing the impact of low quality gene families on the likelihood of each root
   1. Use _DTL_ratio_analysis_ML_diff.py_ create plots sequentially removing the worst quality gene families and resumming the likelihoods for each root, based on either the loss/speciation, duplication/speciation or transfer/speciation rates ratio. 
      - This script requires the following libraries to be installed: sys, glob, os, time, subprocess, pandas, numpy, seaborn and matplotlib
      - The script requires two text files named "roots_to_test.txt", and "species_list.txt". The roots_to_test.txt file should contain the names of the directories of the rec_files. The species_list.txt file should contain all the species in the dataset and should match the formatting of species tree used in the reconciliation analysis.
      - The script also requires two command line arguments: The name of the directory of the maximum-likelihood root determined by an AU test and the metric you want to remove gene families by D,T,L, LS, DS or TS.
      - The script also remakes the plots only using high species representation reconciliation files. The species representation cutoff is currently set to %50 - this can be manually changed by altering the make_high_rep_df() function.
      - Run the script in the same directory as ‘write_consel_file.py’
      - Run the script with the following command to assess the LS ratio under Ecoli as the most likely root: _python3 DTL_ratio_analysis_ML_diff.py Ecoli LS_
      - Plots and data produced will be saved to a new directory e.g. ‘LS_ratio_results’. The directory name will change with the metric used. 
      - Plots will be saved as both png and svg figures. 

1. Gene content evolution on the most likely rooted species tree
   1. Once the most likely root has been identified, this technique allows users to quantify the relative contributions of duplication, transfer, loss and  origination in the gene content evolution.  These predicted values can be calculated for every branch of the species tree, enabling insights into how gene content evolves over time. Here we provide a script that calculates the predicted total number of duplication, transfer, loss and origination events for all branches of the species tree and estimates the number of genes present in the ancestral genome at each internal tree node.
       - Use the script _branchwise_number_of_events.py_ > dtlo.tsv
       - The script requires no command line arguments
       - Run the script in the directory of the root containing the uml_rec you want to analyze.
       - Call the script with _python3 branchwise_number_of_events.py > dtlo.tsv_
       - The output will be a tab separated file. 
 
1. Reconstructing ancestral nodes
   1. The ALE output also provides estimates of the gene families present at each node. Therefore we can model the presence and absence of gene families at internal nodes - reconstructing ancestral genomes. We provide the following script to reconstruct all the genes present at internal nodes. The command line argument is the copy number cutoff - to make the ancestral reconstruction more stringent increase the copy number. 
        - Use the script _Ancestral_reconstruction_copy_number_.py 0.5_ 
        - This script creates separate csv files for each node. 
        - In each of the csv files will be a list of the gene families at that node. 





