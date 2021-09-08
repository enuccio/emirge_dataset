from numpy import average, median
from Bio import SeqIO
import re

### USAGE: Slice a gapped Silva alignment to 1300bp of aligned nucleotides
### Python 2.7


### UPDATE THESE FILENAMES
filename = "combined_seqs_aligned.fasta" #input file to be sliced
out_filename = "combined_seqs_aligned_sliced.fasta"
outfh = open(out_filename,"w")


def seq_length(seq):
    return seq.count("A") + seq.count("T") + seq.count("G") + seq.count("C") + seq.count("U")

def for_rev_gaps(sequence, pattern):
    match_for = re.search(pattern, sequence)
    match_rev = re.search(pattern, sequence[::-1])
    forward_len = len(match_for.group())
    reverse_len = len(match_rev.group())
    return forward_len, reverse_len

def find_gap_lengths(filename, min_bp, max_gap):
    fh = open(filename, "U")
    pattern = re.compile("^[-.]*")
    forward_gaps = []
    reverse_gaps = []
    #length_list = []
    idx = 0
    
    for record in SeqIO.parse(filename, "fasta"):
        header = record.id
        sequence = str(record.seq)
        
        forward_len, reverse_len = for_rev_gaps(sequence, pattern)
        seq_len = seq_length(sequence)
        
        #length_list.append(seq_length)

        if forward_len > max_gap or reverse_len > max_gap:
            #print header, forward_len, reverse_len, seq_length(sequence)
            idx +=1 ## if you want to know how many sequences you're losing
            continue
        else:
            forward_gaps.append(forward_len)
            reverse_gaps.append(reverse_len)
            
    return max(forward_gaps), max(reverse_gaps), median(reverse_gaps)

for gaps in range(3500,3501): # If need to figure out how many max_gaps are allowable
    min_bp = 0 # Change if you want to remove sequences that are shorter than a particular length
    max_gap = gaps

    forward_gaps, reverse_gaps, median_gaps = find_gap_lengths(filename, min_bp, max_gap)
    pattern = re.compile("^[-.]*")

    idx = 0
    lost_gaps = 0
    lost_seqs = 0
    both_lost = 0

    for record in SeqIO.parse(filename, "fasta"):
        header = record.id
        sequence = str(record.seq)

        sliced = sequence[forward_gaps:-reverse_gaps]
        seq_len = seq_length(sequence)
        seq_len_sliced = seq_length(sliced)

        forward_len, reverse_len = for_rev_gaps(sequence, pattern)

        if seq_len_sliced < min_bp:
            #print "Rejected: ", seq_len_sliced, forward_len, reverse_len
            lost_seqs +=1

        if (forward_len > max_gap) or (reverse_len > max_gap):
            lost_gaps +=1
        else:
            #sliced = sequence[6258:-5598]
            #print "Sequence length: ", seq_len_sliced
            idx +=1

        #|perc_sim, seq_id, silva_id, align_length, mismatch, evalue, bitscore, largest_match = blast_dict[header]
        #|taxonomy = taxonomy_dict[silva_id]

        if (forward_len > max_gap) and (seq_len >= min_bp):
            pass
            #print record.id, "Left-hand gaps, but long sequence:", seq_len
            #print "%s percent similar (over %s bases) to %s" % (perc_sim, align_length, taxonomy)


        if (reverse_len > max_gap) and (seq_len >= min_bp):
            #print reverse_len, record.id, "Right-hand gaps, but long sequence:", seq_len, perc_sim, taxonomy
            pass
        
        if (forward_len <= max_gap) and (reverse_len <= max_gap) and (seq_len >= min_bp):
            outfh.write(">"+ header +"\n")
            outfh.write(sliced +"\n")


        '''
        if seq_len > 1500:
            print "forward lost:", seq_length(linestrip[:forward_gaps])
            print "reverse lost:", seq_length(linestrip[-reverse_gaps:])
        '''
        

    print "%i allowed gaps, Forward gaps: %i, Reverse gaps: %i, Median gaps: %i, Lost due to gaps: %i, Lost due to length: %i, Kept: %i" % (gaps, forward_gaps, reverse_gaps, median_gaps, lost_gaps, lost_seqs, idx)


#fh.close()
#outfh.close()
