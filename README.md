# The Distribution of Several Genomic Virulence Determinants Does Not Corroborate the Established Serotyping Classification of *Bacillus thuringiensis*.
A repository with working scripts for IJMS (MDPI) 2020 paper


<img src="https://github.com/lab7arriam/IJMS_2020/blob/main/pics/Fig3.svg?sanitize=true">
<img src="https://github.com/lab7arriam/IJMS_2020/blob/main/pics/Fig5.svg?sanitize=true">

## Reference 

Anton E. Shikov, Yury V. Malovichko  *et al.*, [The Distribution of Several Genomic Virulence Determinants Does Not Corroborate the Established Serotyping Classification of *Bacillus thuringiensis*](https://www.mdpi.com/journal/ijms).
## Contents 

This repository contains scripts used for data preparation for the manuscript. Please consult the Methods section in the paper for extra details. 

## Scripts' description
<ul>
  <li><em>agregate_flagelin_data.py</em>, <em>compare_flagellin_sets.py</em> – aggregating the distribution of lengths and abundances of the flagellin sequences clustered via Roary (Table S10 in the article);</li>
  <li><em>annotate_bt_assemblies.py</em> – summarizing metadata of the assemblies inspected (Table S4 in the article);</li>
  <li><em>bt_pangenome_tree.R</em> – visualizing trees and heatmaps, performing PCA;</li>
  <li><em>calculate_mash_dist.py</em> – constructing a heatmap with paired mash distance scores for all analyzed genomes (Table S12 in the article);</li>
  <li><em>calculate_mean_support.py</em> – calculating mean supporting values for phylogenetic trees in Newick format;</li>
  <li><em>check_lengths.py</em>, <em>CheckLengths.py</em> - calculates sequence length for each sequence in the <i>fasta</i> and prints it in <b>sequence ID - sequence length</b> notation;</li>
  <li><em>check_lengths_massive.py</em> - performs <i>check_lengths.py</i> over a directory containing <i>multifasta</i> files solely;</li>
  <li><em>cluster_pivot_PCR.py</em> - adds initial cluster names and cluster names for the obtained amplicons (not referred to in the final version of the manuscript);</li>
  <li><em>compHmmToAnnot.py</em> - parses a table containing names of Roary-deduced orthologs, then compares them to the entries stored in the HMM output folder and extracts sequences matching the original cluster sequence names from the HMM outputs folder to a separate directory;</li>
  <li><em>compare_genomes_full.py</em> - constructing a heatmap with paired genome identity values using minimap2 for all analyzed genomes (Table S12 in the article);</li>
  <li><em>download_hags.py</em> - downloads the <i>hag</i> gene sequences from Xu and Côté, 2006 (DOI: 10.1128/AEM.00328-06);</li>
  <li><em>downloading.sh</em> – downloading Bt assemblies the NCBI assembly database;</li>
  <li><em>extract_proteins_for_trees.py</em> – extracting protein sequences from Roary-emanated pangenome for a specific gene cluster;</li>
  <li><em>ExtractByClusterName.py</em> - extracts nucleotide sequences by the accessions stored in the Roary cluster table;</li>
  <li><em>ExtractByLength.py</em> - finds the longest/shortest sequence in the sequence file</em></li>
  <li><em>extract_proteins.py'</em> - </li>
  <li><em>get_mean_identity.py</em> – evaluating mean paired sequence identity in the fasta file;</li>
  <li><em>parse_tree_topology.py</em> – assessing the lengths of subtrees containing representatives of <em>Bt</em> serovars (table S14 in the article);</li>
  <li><em>roary_stat.py</em>, <em>summ_roary_stat.py</em> – excluding assemblies from the Roaty-generated pangenome based on the abundance of common gene clusters;</li>
  <li><em>summarize_proteins_from_dige.py</em> – gathering the gene presence/absence results based on diamond blastp results (table S6 in the article).</li>
</ul>
