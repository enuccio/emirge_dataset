import re
import os
from math import ceil

### GOAL: Take trimmed alignment file (from trim_alignment.py) and produce a dereplicated file suitable for UPARSE -cluster_otus
### NOTE: Was written with python 2.7

########## UPDATE THESE FILES NAMES:
# 1. Headers from combined_seqs.fna:
file1 = "qiime_headers_prok.txt"
fh_qiime = open(file1, "U")

# 2. Mapping file from add_qiime_labels.py:
file2 = "id_map_emirge.txt"
fh_mapping = open(file2, "U")

# 3. Bowtie stats file
file3 = "bowtie_prok_stats.txt"
fh_bowtie = open(file3, "U")

# 4. Trimmed alignment file from trim_alignment.py
file4 = "combined_seqs_aligned_sliced_1300bp.fasta"
fh_align = open(file4, "U")

# Output file...use whatever name you want:
file5 = "combined_seqs_aligned_sliced_1300bp_dereplicated.fa"
outfh = open(file5, "w")
##########


dict_qiime = {}
dict_emirge = {}
dict_relabund = {}
dict_taxonomy = {}
dict_bowtie = {}
run_id_list = []
min_seq = 1000 #arbitrary starting number
    
# Parse QIIME headers.  Create dictionary where qiimeID (eg, R1_1000001) maps to relative abundance of sequence (eg, 0.012019)
for line in fh_qiime:
    linesplit = re.split('>|\s|=',line)
    emirgeID = linesplit[2] #didn't use this
    qiimeID = linesplit[1]
    relabund = float(linesplit[8]) # 4 is Prior and 8 is NormPrior...decided on NormPrior
    dict_qiime[qiimeID] = relabund

# Parse bowtie results for number of sequences. Create dictionary where runID (eg, B1) maps to number of seqs per emirge run (eg, 603898)
for line in fh_bowtie:
    if line[0] == "/":
        filename = os.path.basename(line.strip()).split("_")
        runID = filename[0]
    else:
        aligned_seqs = line.strip().split(" ")[0]
        dict_bowtie[runID] = float(aligned_seqs)

# Parse alignment header for qiimeID, calc number of seqs for that qiimeID, put in UPARSE format
# Also remove periods and dashes from alignment file
# Example header from alignment file: >100050 B1_1038239 1..1380
# UPARSE format: >1000741_1;size=6;
pattern = re.compile("[.-]")
for line in fh_align:
    if line[0] == ">":
        linesplit = line.split()
        qiimeID = linesplit[0][1:] # new file, eg: ['>R1_1000000', '3880.R1.148846', 'Prior=0.016292', 'Length=1490', 'NormPrior=0.015459', '1..1490'
        #qiimeID = linesplit[1] # for the original qiime alignment file
        relabund = dict_qiime[qiimeID]
        runID = qiimeID.split("_")[0]
        aligned_seqs = dict_bowtie[runID]
        seq_number = relabund*aligned_seqs
        new_header = ">" + qiimeID +";size="+ str(int(ceil(seq_number))) +";\n"

        if seq_number < min_seq: #Find the minimum number of sequences...will call these "singletons"
            min_seq = seq_number
    else:
        sequence = pattern.sub("",line)
        outfh.write(new_header)
        outfh.write(sequence)
        #outfh.write(line)

print str(int(ceil(min_seq)))
