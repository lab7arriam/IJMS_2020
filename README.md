# The Distribution of Several Genomic Virulence Determinants Does Not Corroborate the Established Serotyping Classification of *Bacillus thuringiensis*.
A repository with working scripts for IJMS (MDPI) 2020 paper

<img src="https://github.com/lab7arriam/IJMS_2020/blob/main/pics/Fig5.svg?sanitize=true">

## Reference 

Anton E. Shikov, Yury V. Malovichko  *et al.*, [The Distribution of Several Genomic Virulence Determinants Does Not Corroborate the Established Serotyping Classification of *Bacillus thuringiensis*](https://www.mdpi.com/journal/ijms).
## Contents 

This repository contains scripts used for data preparation for the manuscript. Please consult the Methods section in the paper for extra details. 

## Scripts' description
<ul>
  <li><em>bt_pangenome_tree.R</em> – visualizing trees and heatmaps, performing PCA;</li>
  <li><em>downloading.sh</em> – downloading Bt assemblies the NCBI assembly database;</li>
  <li><em>annotate_bt_assemblies.py</em> – summarizing metadata of the assemblies inspected (Table S4 in the article);</li>
  <li><em>agregate_flagelin_data.py</em>, <em>compare_flagellin_sets.py</em> – aggregating the distribution of lengths and abundances of the flagellin sequences clustered via Roary (Table S10 in the article);</li>
  <li><em>calculate_mash_dist.py</em> – constructing a heatmap with paired mash distance scores for all analyzed genomes (Table S12 in the article);</li>
  <li><em>calculate_mean_support.py</em> – calculating mean supporting values for phylogenetic trees in Newick format;</li>
  <li><em>compare_genomes_full.py</em>– constructing a heatmap with paired genome identity values using minimap2 for all analyzed genomes (Table S12 in the article);</li>
  <li><em>extract_proteins_for_trees.py</em> – extracting protein sequences from Roary-emanated pangenome for a specific gene cluster;</li>
  <li><em>get_mean_identity.py</em> – evaluating mean paired sequence identity in the fasta file;</li>
  <li><em>parse_tree_topology.py</em> – assessing the lengths of subtrees containing representatives of <em>Bt</em> serovars (table S14 in the article);</li>
  <li><em>roary_stat.py</em>, <em>summ_roary_stat.py</em> – excluding assemblies from the Roaty-generated pangenome based on the abundance of common gene clusters;</li>
  <li><em>summarize_proteins_from_dige.py</em> – gathering the gene presence/absence results based on diamond blastp results (table S6 in the article).</li>
</ul>
