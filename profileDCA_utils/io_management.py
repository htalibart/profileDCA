import shutil, pathlib, json, csv
import numpy as np
import pandas as pd
from Bio import AlignIO, SeqIO

from profileDCA_utils.global_variables import ALPHABET, ALPHABET_DICT
from profileDCA_utils import ppfunctions

def get_int_msa_array(msa_file, alphabet_dict=ALPHABET_DICT):
    """ returns MSA in (N,L) np.array format where MSA[n,i] = state of letter in sequence n in column in in multiple sequence alignment in msa_file"""
    msa = AlignIO.read(str(msa_file), 'fasta')
    L = msa.get_alignment_length()
    N = msa.__len__()
    int_array = np.zeros((N,L), dtype=int)
    for i in range(L):
        for n in range(N):
            int_array[n,i] = ppfunctions.get_letter_index(msa[n,i], alphabet_dict)
    return int_array

def mrf_to_files(mrf, output_folder):
    for key in mrf:
        np.save(str(output_folder/key), mrf[key])

def mrf_from_folder(potts_folder):
    mrf = {}
    for key in ['v','w','v_full', 'mrf_pos_to_seq_pos']:
        mrf[key] = np.load(str(potts_folder/(key+'.npy')))
    return mrf

def get_first_sequence_in_fasta_file(seq_file):
    return str(list(SeqIO.parse(str(seq_file), "fasta"))[0].seq)

def get_first_sequence_name(seq_file):
    return str(list(SeqIO.parse(str(seq_file), "fasta"))[0].id)

def get_first_record_sequence_and_name(seq_file):
    first_record = list(SeqIO.parse(str(seq_file), "fasta"))[0]
    return str(first_record.seq), str(first_record.id)


def copy_file(origin, destination):
    if origin!=destination:
        shutil.copy(str(origin), str(destination))

def write_readme(folder, **kwargs):
    p = folder/'README.txt'
    with p.open(mode='w') as f:
        json.dump(kwargs, f, default=str)

def get_parameters_from_readme_file(readme_file):
    params = json.load(open(str(readme_file)))
    return params

def dict_to_csv(d, input_csv_file):
    with open(str(input_csv_file), 'w') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(list(d.keys()))
        csv_writer.writerow(list(d.values()))

def write_positions_to_csv(positions_list, output_file):
    """ PPalign list of aligned positions to csv file """
    with open(str(output_file), 'w') as f:
        csvwriter = csv.writer(f)
        csvwriter.writerow(['pos_ref','pos_2'])
        for ind in range(len(positions_list[0])):
            row = [positions_list[seq_index][ind] for seq_index in range(2)]
            csvwriter.writerow(row)

def get_infos_solver_dict_from_ppalign_output_file(infos_res_file):
    df = pd.read_csv(infos_res_file)
    return df.loc[0].to_dict()

def get_aligned_positions_from_ppalign_output_file(aln_res_file):
    """ get {"pos_ref":list of aligned positions in first Potts model, "pos_2":  list of aligned positions in second Potts model} from PPalign output .csv file """
    with open(str(aln_res_file), 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        header = next(csv_reader)
        aln_dict = {name:[] for name in header}
        for row in csv_reader:
            for k in range(2):
                if (row[k]=='-'):
                    value='-'
                else:
                    try:
                        value = int(row[k])
                    except Exception as e:
                        value = None
                aln_dict[header[k]].append(value)
    assert(len(aln_dict['pos_ref'])==len(aln_dict['pos_2']))
    return [aln_dict['pos_ref'], aln_dict['pos_2']]


def get_aligned_positions_with_gaps_from_ppalign_output_file(aln_with_gaps_res_file):
    """ get {"pos_ref":list of aligned positions with gaps in first Potts model, "pos_2":  list of aligned positions with gaps in second Potts model} from PPalign output .csv file """
    with open(str(aln_with_gaps_res_file), 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        header = next(csv_reader)
        aln_with_gaps_dict = {name:[] for name in header}
        for row in csv_reader:
            for k in range(2):
                if (row[k]=='-'):
                    value='-'
                else:
                    try:
                        value = int(row[k])
                    except Exception as e:
                        value = None
                aln_with_gaps_dict[header[k]].append(value)
    assert(len(aln_with_gaps_dict['pos_ref'])==len(aln_with_gaps_dict['pos_2']))
    return [aln_with_gaps_dict['pos_ref'], aln_with_gaps_dict['pos_2']]

def create_folder(folder_path):
    if not folder_path.is_dir():
        folder_path.mkdir()
    else:
        raise Exception("Folder already exists")

