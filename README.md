# emirge_dataset
Scripts used in Nuccio et al. 2021 manuscript. Note that the scripts are provided for archival purposes and not maintained.

## sequence_analysis.sh
The pipeline script "sequence_analysis.sh" provides code used to create EMIRGE assembled sequences, as well as run QIIME and UPARSE. At Step 4, the script uses two custom python scripts created to trim the aligned EMIRGE sequences containing gaps to 1300 nucleotides (trim_alignment.py), as well as convert the trimmed alignments into a format readable by UPARSE that contains read abundances calculated from EMIRGE and bowtie (emirge_relabund_to_uparse.py).
