import numpy as np

from Bio import SeqIO

from profileDCA_utils import io_management as iom
from profileDCA_utils import global_variables, ppfunctions

import pathlib, tempfile, os, subprocess, shutil

def copy_file(origin, destination):
    shutil.copy(str(origin), str(destination))

def clean_fasta_file(input_fasta, output_fasta, alphabet=global_variables.ALPHABET):
    # remove sequences with letters not in alphabet and capitalize everything
    all_records = list(SeqIO.parse(str(input_fasta), 'fasta'))
    clean_records = []
    for record in all_records:
        if all(letter in alphabet for letter in str(record.seq).upper()):
            clean_records.append(record.upper())
    SeqIO.write(clean_records, str(output_fasta), 'fasta')

def call_hhfilter(input_file, output_file, hhid):
    cmd = "hhfilter -i "+str(input_file)+" -o "+str(output_file)+" -id "+str(hhid)
    subprocess.Popen(cmd, shell=True).wait()
    if not output_file.is_file():
        raise Exception("HHfilter failed")

def make_fasta_nonredundant(input_fasta, output_fasta, seq_id_threshold):
    hhfilter_threshold = int(100*seq_id_threshold)
    call_hhfilter(input_fasta, output_fasta, hhfilter_threshold)


def create_file_with_less_sequences(input_file, output_file, nb_sequences, fileformat="fasta"):
    records = list(SeqIO.parse(str(input_file),fileformat))
    with open(str(output_file), 'w') as f:
        SeqIO.write(records[:nb_sequences], f, fileformat)

def process_msa_for_inference(input_msa_file, output_msa_file, seq_id_threshold=0.8, max_nb_sequences=None):
    iom.copy_file(input_msa_file, output_msa_file)
    clean_fasta_file(input_msa_file, output_msa_file)
    if seq_id_threshold<1:
        make_fasta_nonredundant(output_msa_file, output_msa_file, seq_id_threshold)
    if max_nb_sequences!=None:
        create_file_with_less_sequences(output_msa_file, output_msa_file, max_nb_sequences)

def build_int_msa_array_for_all_sequence_positions(msa_file, sequence_file, alphabet_dict=global_variables.ALPHABET_DICT):
    original_msa_int_array = iom.get_int_msa_array(msa_file)
    N, L_msa = original_msa_int_array.shape
    sequence = iom.get_first_sequence_in_fasta_file(sequence_file)
    L_seq = len(sequence)
    if L_msa==L_seq:
        return original_msa_int_array
    else:
        seq_int_array = np.zeros((N,L_seq))
        msa_first_seq = iom.get_first_sequence_in_fasta_file(msa_file)
        seq_pos_to_msa_pos = ppfunctions.get_pos_first_seq_to_second_seq(sequence, msa_first_seq)
        for seq_pos in range(L_seq):
            msa_pos = seq_pos_to_msa_pos[seq_pos]
            if (msa_pos!=None):
                seq_int_array[:,seq_pos] = original_msa_int_array[:,msa_pos]
            else:
                seq_int_array[:,seq_pos] = global_variables.ALPHABET_DICT['-']
                seq_int_array[0,seq_pos] = pputils.get_letter_index(self.sequence[seq_pos], alphabet_dict)
        return seq_int_array
