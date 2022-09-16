from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def get_aligned_sequence_positions_with_gaps(aligned_mrf_positions_with_gaps, mrf_pos_to_seq_pos_list):
    aligned_sequence_positions_with_gaps = [[],[]]
    aln_length = len(aligned_mrf_positions_with_gaps[0])
    seq_pos_previously_aligned = [-1,-1]
    for pos_in_aln in range(aln_length):
        current_seq_pos = []
        for seq_index in range(2):
            other_seq_index = (seq_index+1) % 2 # index of other sequence in list
            if aligned_mrf_positions_with_gaps[seq_index][pos_in_aln]=='-':
                seq_pos = '-'
            else:
                seq_pos = mrf_pos_to_seq_pos_list[seq_index][aligned_mrf_positions_with_gaps[seq_index][pos_in_aln]]
                if seq_pos!=None:
                    for previous_pos in range(seq_pos_previously_aligned[seq_index]+1, seq_pos):
                        aligned_sequence_positions_with_gaps[seq_index].append(previous_pos)
                        aligned_sequence_positions_with_gaps[other_seq_index].append('-')
                    seq_pos_previously_aligned[seq_index] = seq_pos
                else:
                    seq_pos='-'
            current_seq_pos.append(seq_pos)
        for seq_index in range(2):
            aligned_sequence_positions_with_gaps[seq_index].append(current_seq_pos[seq_index])
    for seq_index in range(2):
        other_seq_index = (seq_index+1) % 2
        last_existing_pos = [pos for pos in mrf_pos_to_seq_pos_list[seq_index] if pos!=None][-1]
        for previous_pos in range(seq_pos_previously_aligned[seq_index]+1, last_existing_pos+1):
            aligned_sequence_positions_with_gaps[seq_index].append(previous_pos)
            aligned_sequence_positions_with_gaps[other_seq_index].append('-')
    return aligned_sequence_positions_with_gaps



def aligned_positions_with_gaps_to_aligned_sequences(aligned_positions_with_gaps, sequences):
    aligned_sequences = ["",""]
    for seq_index in range(2):
        for pos_in_seq in aligned_positions_with_gaps[seq_index]:
            if pos_in_seq == '-':
                aligned_sequences[seq_index] += '-'
            else:
                aligned_sequences[seq_index] += sequences[seq_index][pos_in_seq]
    return aligned_sequences


def aligned_positions_with_gaps_to_aligned_sequences_with_lowercase(all_aligned_positions_with_gaps, confident_aligned_positions, sequences):
    aligned_sequences = ["",""]
    for seq_index in range(2):
        index_in_confident=0
        for pos_in_seq in all_aligned_positions_with_gaps[seq_index]:
            if pos_in_seq == '-':
                aligned_sequences[seq_index] += '-'
            else:
                if (index_in_confident<len(confident_aligned_positions[seq_index])):
                    if pos_in_seq==confident_aligned_positions[seq_index][index_in_confident]:
                        aligned_sequences[seq_index] += sequences[seq_index][pos_in_seq].upper()
                        index_in_confident+=1
                    else:
                        aligned_sequences[seq_index] += sequences[seq_index][pos_in_seq].lower()
                else:
                    aligned_sequences[seq_index] += sequences[seq_index][pos_in_seq].lower()
    return aligned_sequences


def translate_aligned_positions_without_gaps(aligned_positions_without_gaps_1, list_of_pos1_to_pos2):
    """
        translates each aligned position in aligned_positions_without_gaps_1 to aligned positions using pos list correspondences list_of_pos1_to_pos2
        resulting alignment is the same length
        it does NOT take gaps into account
    """
    aligned_positions_without_gaps_2 = [[],[]]
    for seq_index in range(2):
        for pos_in_aln in range(len(aligned_positions_without_gaps_1[0])):
            pos_1 = aligned_positions_without_gaps_1[seq_index][pos_in_aln]
            pos_2 = list_of_pos1_to_pos2[seq_index][pos_1]
            aligned_positions_without_gaps_2[seq_index].append(pos_2)
    return aligned_positions_without_gaps_2



def write_aligned_sequences_in_fasta_file(aligned_sequence_positions_with_gaps, sequences, names, output_fasta_file):
    seqs_aligned = aligned_positions_with_gaps_to_aligned_sequences(aligned_sequence_positions_with_gaps, sequences)
    seq_records = [SeqRecord(Seq(s), id=name, description='') for s,name in zip(seqs_aligned, names)]
    with open(str(output_fasta_file), 'w') as f:
        SeqIO.write(seq_records, f, "fasta")

def write_aligned_sequences_in_fasta_file_with_lowercase(aligned_sequence_positions_with_gaps, confident_aligned_sequence_positions, sequences, names, output_fasta_file):
    seqs_aligned = aligned_positions_with_gaps_to_aligned_sequences_with_lowercase(aligned_sequence_positions_with_gaps, confident_aligned_sequence_positions, sequences)
    seq_records = [SeqRecord(Seq(s), id=name, description='') for s,name in zip(seqs_aligned, names)]
    with open(str(output_fasta_file), 'w') as f:
        SeqIO.write(seq_records, f, "fasta")

def reverse_pos1_to_pos2(pos1_to_pos2, L2):
    """
        pos2_to_pos1[pos2] = index of pos2 in pos1_to_pos2 if it exists, else None
        needs length @L2 of second model to complete resulting list with remaining positions
    """
    pos2_to_pos1 = []
    previous_pos2 = -1
    for pos1 in range(len(pos1_to_pos2)):
        current_pos2 = pos1_to_pos2[pos1]
        if current_pos2!=None:
            for none_pos in range(previous_pos2+1,current_pos2):
                pos2_to_pos1.append(None)
            previous_pos2 = current_pos2
            pos2_to_pos1.append(pos1)
    for pos1 in range(len(pos2_to_pos1),L2):
        pos2_to_pos1.append(None)
    return pos2_to_pos1

def get_insert_opens_for_mrf_pos_to_seq_pos(gap_open, mrf_pos_to_seq_pos, no_free_end_gaps=False):
    """ returns a list of length self.potts_model.ncol+1 of open insertion penalties where penalty is set to 0 if positions where trimmed after and @gap_open otherwise """
    if no_free_end_gaps:
        insert_opens = [gap_open]
    else:
        insert_opens = [0]
    previous_seq_pos = 0
    for mrf_pos in range(1,len(mrf_pos_to_seq_pos)):
        seq_pos = mrf_pos_to_seq_pos[mrf_pos]
        if (seq_pos==None):
            insert_opens.append(gap_open)
        else:
            if seq_pos-previous_seq_pos>1:
                insert_opens.append(0)
            else:
                insert_opens.append(gap_open)
            previous_seq_pos = seq_pos
    if no_free_end_gaps:
        insert_opens.append(gap_open)
    else:
        insert_opens.append(0)
    return insert_opens
