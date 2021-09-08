##### Sequence analysis for EMIRGE dataset #####
# Software versions: USEARCH v7, QIIME 1.8, python 2.7


### STEP 1: Generate EMIRGE sequences (run for each sample)
# Run trimmomatic on paired reads
java -classpath /opt/trimmomatic-0.22.jar org.usadellab.trimmomatic.TrimmomaticPE -phred33 sample1_R1_001.fastq sample2_R2_001.fastq R1_R1_paired_Q3_C10.fastq R1_R1_unpaired_Q3_C10.fastq R1_R2_paired_Q3_C10.fastq R1_R2_unpaired_Q3_C10.fastq TRAILING:3 MINLEN:60 HEADCROP:10

# Run emirgy_amplicon.py (example)
nice /usr/local/bin/emirge_amplicon.py /directory_with_sequences/sample1 -1 R1_R1_paired_Q3_C10.fastq -2 R1_R1_unpaired_Q3_C10.fastq -i 234 -s 85 -f SSU_Silva_Euk_97.fasta -b SSU_Silva_Euk_97_fix_db_btindex -l 151 -a12 --phred33

# Run emirge_rename_fasta.py - rewrites an emirge fasta file to include proper sequence names and prior probabilities (abundance estimates)
# Note: find in folder containing the different iterations
emirge_rename_fasta.py iter.40 > sample1_iter.40.fasta

# Put sequences in QIIME format
# Note: Rename OTU IDs in EMIRGE files so they're unique (this is important for creating the OTU table later)
# Replace each | with appropriate run name (eg, .B1.)
sed 's/|/.B1euk./' sample1_iter.40.fasta > sample1_unique_iter.40.fasta

# Create a mapping file
euk_id_map_emirge.txt # make a text file that assigns a sample ID to a fasta file, for example:
'''
sample1	sample1_unique_iter.40.fasta
sample2	sample2_unique_iter.40.fasta
sample2	sample3_unique_iter.40.fasta
# etc
'''

# Add QIIME labels and create combined_seqs.fna using add_qiime_labels.py
# Note: This QIIME script takes a directory, a metadata mapping file, and a column name that contains the fasta file names that SampleIDs are associated with, combines all files that have valid fasta extensions into a single fasta file, with valid QIIME fasta labels
add_qiime_labels.py -i /path_to_sequence_directory/ -m euk_id_map_emirge.txt -n 1000000 -o /path_to_sequence_directory/


### STEP 2: Set bash variables for sequence analysis
# Path to USEARCH binary
u=/path_to/usearch_v7

# Path to working working directory
d=/path_to/working_directory


### STEP 3: QIIME Align sequences for prokaryotes and eukaryotes
# Prokaryotes (Reference: QIIME 1.8 core_set_aligned.fasta.imputed)
nice align_seqs.py -i $d/combined_seqs.fna -o $d/combined_seqs_aligned.fna --template_fp /opt/qiime1.8/lib/qiime_test_data/align_seqs/core_set_aligned.fasta.imputed --alignment_method pynast --pairwise_alignment_method uclust --min_percent_id 75.0 --min_length 75 

# Eukaryotes (Reference: Silva 111 clustered at 97% similarity, eukaryotes only)
nice align_seqs.py -i $d/combined_seqs.fna -o $d/combined_seqs_aligned.fna --template_fp $d/97_Silva_111_rep_set_euk_aligned_compressed.fasta --alignment_method pynast --pairwise_alignment_method uclust --min_percent_id 75.0 --min_length 75 


### STEP 4: Run python files to prepare alignments for UPARSE
# Trims gapped aligments to 1300bp
python trim_alignment.py

# Produces deprelication file for UPARSE with the number of sequences calulated from bowtie results
python emirge_relabund_to_uparse.py


### STEP 5: Run UPARSE (script hereafter is for prokaryotes; notes provided where relevant for eukaryotes)
# Abundance sort and discard singletons
# Note: Created own dereplication file
# Note: Sequences were calculated.  Minimum number of "seqs" is 3, so 3 == "singleton"
# Note: GAPS NOT ALLOWED

$u -sortbysize $d/combined_seqs_rep_set_aligned_sliced_1300bp_dereplicated.fa -output $d/combined_seqs_rep_set_aligned_sliced_1300bp_sorted.fa -minsize 4

# OTU clustering
$u -cluster_otus $d/combined_seqs_rep_set_aligned_sliced_1300bp_sorted.fa -otus $d/combined_seqs_rep_set_aligned_sliced_1300bp_otus.fa

# Chimera filtering using reference database
# Note: Used gg_99_13_5_centroids.fasta for db, which just contains the sequences of the centroids (about 1/6th of size, used ~2.4 of 4 GB)
# Eukaryotes: Used Silva 119 at clustered at 97% to make reference file smaller
$u -uchime_ref $d/combined_seqs_aligned_sliced_1300bp_otus.fa -db $d/gg_99_13_5_centroids.fasta  -strand plus -nonchimeras $d/combined_seqs_aligned_sliced_1300bp_otus_no_chimeras.fa

# Label OTU sequences OTU_1, OTU_2…
# Note: Update the path to the usearch folder for fasta_number.py
python /path_to_usearch_folder/fasta_number.py $d/combined_seqs_aligned_sliced_1300bp_otus_no_chimeras.fa OTU_ > $d/combined_seqs_aligned_sliced_1300bp_otus_no_chimeras_otulabels.fa

# Map reads (including singletons) back to OTUs
$u -usearch_global $d/combined_seqs_aligned_sliced_1300bp_dereplicated.fa -db $d/combined_seqs_aligned_sliced_1300bp_otus_no_chimeras_otulabels.fa -strand both -id 0.97 -maxaccepts 2000 -uc $d/combined_seqs_aligned_sliced_1300bp_otus_no_chimeras_map_allhits.uc -uc_allhits

# Convert OTU table to biom table
biom convert -i $d/otu_table_emirge_prok_uparse_1300bp.txt -o $d/otu_table_emirge_prok_uparse_1300bp.biom --table-type "otu table"

# Assign taxonomy
# Note: Added header label #OTUID	taxonomy	confidence
# Eukaryotes: reference_seqs_fp is SILVA_119_SSURef_Nr99_tax_silva_centroids.fasta
nice assign_taxonomy.py -o $d/rdp_assigned_taxonomy -i $d/combined_seqs_aligned_sliced_1300bp_otus_no_chimeras_otulabels.fa --reference_seqs_fp $d/gg_13_5.fasta --id_to_taxonomy_fp $d/gg_13_5_taxonomy.txt --assignment_method rdp --confidence 0.8 --rdp_max_memory 28672

# Add taxa to biom table
biom add-metadata -i $d/otu_table_emirge_prok_uparse_1300bp.biom -o $d/otu_table_emirge_prok_uparse_1300bp_w_tax.biom --observation-metadata-fp $d/rdp_assigned_taxonomy/combined_seqs_aligned_sliced_1300bp_otus_no_chimeras_otulabels_tax_assignments.txt

# Remove OTUs only in one sample
filter_otus_from_otu_table.py -i $d/otu_table_emirge_prok_uparse_1300bp_w_tax.biom -o $d/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered.biom -n 2 -s 2

# Summarize biom stats, determine number of sequences for rarification
biom summarize-table -i  $d/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered.biom -o $d/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered.stats


### STEP 6: reate tree using FastTree: Use a 90% similarly tree to constrain tree (constraining tree with a reference tree helps with placement of branches; this is useful for short sequences but also works with long sequences)
# Make fasta file for rep set 
filter_fasta.py -f $d/combined_seqs_aligned_sliced_1300bp_otus_no_chimeras_otulabels.fa -o $d/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered.fasta -b $d/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered.biom

# Cluster OTUs at 90%
nice $u -cluster_fast $d/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered.fasta -id 0.9 -centroids $d/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered_centroids.fasta

# Make alignment for the entire rep set
nice align_seqs.py -i $d/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered.fasta  -o $d/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered_aligned  --template_fp /opt/qiime1.8/lib/qiime_test_data/align_seqs/core_set_aligned.fasta.imputed --alignment_method pynast --pairwise_alignment_method uclust --min_percent_id 75.0 --min_length 75 

# Filter alignment using lanemask
filter_alignment.py -i $d/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered_aligned/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered_aligned.fasta -m /opt/qiime1.8/lib/qiime_test_data/filter_alignment/lanemask_in_1s_and_0s -o $d/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered_aligned/

# Make fasta file of pfiltered aligned 90% clusters
filter_fasta.py -f $d/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered_aligned/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered_aligned_pfiltered.fasta -o $d/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered_aligned/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered_aligned_pfiltered_90perc.fasta -a $d/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered_centroids.fasta

# Make 90% tree (constraint tree)
FastTree -nt -gtr -gamma $d/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered_aligned/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered_aligned_pfiltered_90perc.fasta > $d/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered_aligned/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered_aligned_pfiltered_90perc.tree

# Make constraints file
perl TreeToConstraints.pl < $d/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered_aligned/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered_aligned_pfiltered_90perc.tree > $d/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered_aligned/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered_aligned_pfiltered_90perc.constraint

# Make rep set tree constrained by 90% tree
FastTree -nt -gtr -gamma -constraints $d/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered_aligned/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered_aligned_pfiltered_90perc.constraint < $d/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered_aligned/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered_aligned_pfiltered.fasta > $d/emirge_prok_uparse_rep_set_90perc_constrained.tree


### STEP 7: Run beta-diversity analysis with rarification
# NEEDED TO ADD SINGLE QUOTES AROUND OTU LABELS, e.g., re.sub(r"(OTU_\d+)", r"'\1'", tree)
beta_diversity_through_plots.py -i $d/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered.biom -o $d/beta_diversity_even121737 -m $d/emirge_mapping.txt -t $d/emirge_prok_uparse_rep_set_90perc_constrained_EDITED.tree -e 121737 -f


### STEP 8: Statistical analysis of dataset
# Rhizosphere-control vs. Rhizosphere-litter: R_RL
group_significance.py -i $d/beta_diversity_even121737/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered_even121737.biom -m $d/emirge_mapping.txt -c R_RL -s parametric_t_test -o $d/ttest_R_RL.txt

# Bulk-control vs. Bulk-litter: B_BL
group_significance.py -i $d/beta_diversity_even121737/otu_table_emirge_prok_uparse_1300bp_w_tax_filtered_even121737.biom -m $d/emirge_mapping.txt -c B_BL -s parametric_t_test -o $d/ttest_B_BL.txt

