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
  <li><em>agregate_cry_data.py</em> - summarizes the results about 3-D cry proteins spectra in the assemblies based on CryProcessor results (<b>Table S15</b> in the article);</li>
  <li><em>agregate_flagelin_data.py</em>, <em>compare_flagellin_sets.py</em> – aggregating the distribution of lengths and abundances of the flagellin sequences clustered via Roary (<b>Table S10</b> in the article);</li>
  <li><em>annotate_bt_assemblies.py</em> – summarizing metadata of the assemblies inspected (<b>Table S4</b> in the article);</li>
  <li><em>bt_pangenome_tree.R</em> – visualizing trees and heatmaps, performing PCA (<b>Figures 3-5</b> in the article);</li>
  <li><em>calculate_mash_dist.py</em> – constructing a heatmap with paired mash distance scores for all analyzed genomes (<b>Table S12</b> and <b>Figure S6a</b> in the article);</li>
  <li><em>calculate_mean_support.py</em> – calculating mean supporting values for phylogenetic trees in Newick format (for <b>Table S8</b> in the article);</li>
  <li><em>check_lengths.py</em>, <em>CheckLengths.py</em> - calculates sequence length for each sequence in the <i>fasta</i> and prints it in <b>sequence ID - sequence length</b> notation;</li>
  <li><em>check_lengths_massive.py</em> - performs <i>check_lengths.py</i> over a directory containing <i>multifasta</i> files solely;</li>
  <li><em>cluster_pivot_PCR.py</em> - adds initial cluster names and cluster names for the obtained amplicons (not referred to in the final version of the manuscript);</li>
  <li><em>compHmmToAnnot.py</em> - parses a table containing names of Roary-deduced orthologs, then compares them to the entries stored in the HMM output folder and extracts sequences matching the original cluster sequence names from the HMM outputs folder to a separate directory;</li>
  <li><em>compare_genomes_full.py</em> - constructing a heatmap with paired genome identity values using minimap2 for all analyzed genomes (<b>Table S12</b> and <b>Figure S6b</b> in the article);</li>
  <li><em>download_hags.py</em> - downloads the <i>hag</i> gene sequences from Xu and Côté, 2006 (DOI: 10.1128/AEM.00328-06);</li>
  <li><em>downloading.sh</em> – downloading <em>Bt</em> assemblies the NCBI assembly database;</li>
  <li><em>extract_proteins_for_trees.py</em> – extracting protein sequences from Roary-emanated pangenome for a specific gene cluster;</li>
  <li><em>extract_proteins.py</em> - extracts sequences from the <i>fasta</i> files based on file containing a list of identifiers</li>
  <li><em>extract_proteins_blast.py: extracts query identifiers from the previously filtered BLAST outputs and uses them to fetch protein sequences from the <i?multifasta</i> files;</li>
  <li><em>ExtractByClusterName.py</em> - extracts nucleotide sequences by the accessions stored in the Roary cluster table and fetches sequences from the cluster <i>fasta</i> files;</li>
  <li><em>ExtractByLength.py</em> - finds the longest/shortest sequence in the sequence file</em></li>
  <li><em>fetch_nucleotide.py'</em>: assigns sequences from the HMMer or BLAST output to the Roary cluster reprsentatives adn</li>
  <li><em>get_mean_identity.py</em> – evaluating mean paired sequence identity in the fasta file (for <b>Table S8</b> in the article);</li>
  <li><em>parse_tree_topology.py</em> – assessing the lengths of subtrees containing representatives of <em>Bt</em> serovars (<b>table S14</b> in the article);</li>
  <li><em>roary_stat.py</em>, <em>summ_roary_stat.py</em> – excluding assemblies from the Roary-generated pangenome based on the abundance of common gene clusters;</li>
  <li><em>summarize_proteins_from_dige.py</em> – gathering the gene presence/absence results based on diamond blastp results (<b>table S6</b> in the article).</li>
</ul>
