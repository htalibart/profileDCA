import numpy as np

import pathlib, tempfile, argparse, sys
from profileDCA_build import msa_processing, mrf_inference, sequence_to_msa
from profileDCA_utils import global_variables
from profileDCA_utils import io_management as iom

def get_sub_v(mrf, pos_begin, pos_end):
    """ returns MRF with v_full from position @pos_begin to position @pos_end included and w=0 """
    q = mrf['v'].shape[1]
    new_L = pos_end-pos_begin+1
    if pos_begin!=None and pos_end!=None:
        if pos_end<pos_begin:
            small_v = np.empty(shape=(0,q))
            small_w = np.empty(shape=(0,0,q,q))
        else:
            small_v = mrf['v_full'][pos_begin:pos_end+1]
            small_w = np.zeros((new_L,new_L,q,q))
    else:
        small_v = np.empty(shape=(0,q))
        small_w = np.empty(shape=(0,0,q,q))            
    return {'v':small_v, 'w': small_w}

def get_mrf_from_processed_msa(msa_file, sequence_file=None, max_gap_v=1, max_gap_w=1, uniform_pc_rate=0.5, apply_covariance_matrix_threshold=True, shrinkage_coeff=0.7, reg_lambda_w=1e-4, alphabet_dict=global_variables.ALPHABET_DICT):
    if sequence_file!=None:
        msa_array = msa_processing.build_int_msa_array_for_all_sequence_positions(msa_file, sequence_file, alphabet_dict=alphabet_dict)
    else:
        msa_array = iom.get_int_msa_array(msa_file)
    mrf = mrf_inference.msa_array_to_mrf(msa_array, max_gap_v=max_gap_v, max_gap_w=max_gap_w, uniform_pc_rate=uniform_pc_rate, apply_covariance_matrix_threshold=apply_covariance_matrix_threshold, shrinkage_coeff=shrinkage_coeff, reg_lambda_w=reg_lambda_w, alphabet_dict=alphabet_dict)
    return mrf

def get_mrf_from_unprocessed_msa(msa_file, output_msa_file, sequence_file=None, max_gap_v=1, max_gap_w=1, uniform_pc_rate=0.5, apply_covariance_matrix_threshold=True, shrinkage_coeff=0.7, reg_lambda_w=1e-4, alphabet_dict=global_variables.ALPHABET_DICT, seq_id_threshold=0.8, max_nb_sequences=None):
    msa_processing.process_msa_for_inference(msa_file, output_msa_file, seq_id_threshold=seq_id_threshold, max_nb_sequences=max_nb_sequences)
    mrf = get_mrf_from_processed_msa(output_msa_file, max_gap_v=max_gap_v, max_gap_w=max_gap_w, uniform_pc_rate=uniform_pc_rate, apply_covariance_matrix_threshold=apply_covariance_matrix_threshold, shrinkage_coeff=shrinkage_coeff, reg_lambda_w=reg_lambda_w, alphabet_dict=alphabet_dict, sequence_file=sequence_file)
    return mrf 


def get_mrf_from_sequence(sequence_file, database, output_msa=None, max_gap_v=1, max_gap_w=1, uniform_pc_rate=0.5, apply_covariance_matrix_threshold=True, shrinkage_coeff=0.7, reg_lambda_w=1e-4, alphabet_dict=global_variables.ALPHABET_DICT, seq_id_threshold=0.8, max_nb_sequences=None):
    if output_msa!=None:
        processed_msa_file = output_msa
    else:
        processed_msa_file = pathlib.Path("/tmp/")/(next(tempfile._get_candidate_names())+'.fasta')
    sequence_to_msa.sequence_to_processed_msa(sequence_file, database, processed_msa_file, seq_id_threshold=seq_id_threshold, max_nb_sequences=max_nb_sequences)
    mrf = get_mrf_from_processed_msa(processed_msa_file, max_gap_v=max_gap_v, max_gap_w=max_gap_w, uniform_pc_rate=uniform_pc_rate, apply_covariance_matrix_threshold=apply_covariance_matrix_threshold, shrinkage_coeff=shrinkage_coeff, reg_lambda_w=reg_lambda_w, alphabet_dict=alphabet_dict, sequence_file=sequence_file)
    if output_msa==None:
        processed_msa_file.unlink()
    return mrf


def build_mrf_to_folder(potts_folder, sequence_file, msa_file=None, database=None, max_gap_v=1, max_gap_w=1, uniform_pc_rate=0.5, apply_covariance_matrix_threshold=True, shrinkage_coeff=0.7, reg_lambda_w=1e-4, alphabet_dict=global_variables.ALPHABET_DICT, seq_id_threshold=0.8, max_nb_sequences=None):
    if not potts_folder.is_dir():
        potts_folder.mkdir()
    iom.copy_file(sequence_file, potts_folder/"sequence.fasta")
    if msa_file==None:
        if database==None:
            raise Exception("Missing database path for HHblits")
        else:
            msa_file = potts_folder/"msa.fasta"
            mrf = get_mrf_from_sequence(sequence_file, database, msa_file, max_gap_v=max_gap_v, max_gap_w=max_gap_w, uniform_pc_rate=uniform_pc_rate, apply_covariance_matrix_threshold=apply_covariance_matrix_threshold, shrinkage_coeff=shrinkage_coeff, reg_lambda_w=reg_lambda_w, alphabet_dict=alphabet_dict, seq_id_threshold=seq_id_threshold, max_nb_sequences=max_nb_sequences)
    else:
        mrf = get_mrf_from_unprocessed_msa(msa_file, potts_folder/"msa.fasta", sequence_file=sequence_file, max_gap_v=max_gap_v, max_gap_w=max_gap_w, uniform_pc_rate=uniform_pc_rate, apply_covariance_matrix_threshold=apply_covariance_matrix_threshold, shrinkage_coeff=shrinkage_coeff, reg_lambda_w=reg_lambda_w, alphabet_dict=alphabet_dict, seq_id_threshold=seq_id_threshold, max_nb_sequences=max_nb_sequences)
    iom.mrf_to_files(mrf, potts_folder)
    return mrf



def main(args=sys.argv[1:]):
    parser = argparse.ArgumentParser()

    # files
    file_args = parser.add_argument_group('file_args')
    file_args.add_argument('-pf', '--potts_folder', help="Output feature folder", type=pathlib.Path)
    file_args.add_argument('-s', '--sequence_file', help="Sequence file", type=pathlib.Path)
    file_args.add_argument('-msa', '--msa_file', help="Alignment file", type=pathlib.Path, default=None)

    # hhblits
    hhblits_args = parser.add_argument_group('hhblits_args')
    hhblits_args.add_argument('-hd', '--database', help="Database path for HHblits call", default=None)
    
    # aln files processing
    aln_processing_args = parser.add_argument_group('aln_processing_args')
    aln_processing_args.add_argument('--seq_id_threshold', help="Sequence identity threshold for HHfilter (between 0 and 1, default:0.8)", type=float, default=0.8)
    aln_processing_args.add_argument('-maxnb', '--max_nb_sequences', help="Max. nb sequences in the alignment (if alignment has more sequences that @max_nb_sequences after filtering and before trimming, all sequences after n° @max_nb_sequences will be deleted from the alignment. Default: 2000)", type=int, default=2000)

    # gap processing
    gap_processing_args = parser.add_argument_group('gap_processing_args')
    gap_processing_args.add_argument('--max_gap_w', help="Columns with more than max_gap_w fraction of gaps will not be considered in coupling inference, their wij will be set to 0 and their vi will be computed as if they were fields of an independent model", type=float, default=0.7)
    gap_processing_args.add_argument('--max_gap_v', help="Columns with more than max_gap_v fraction of gaps will not be in the model", type=float, default=0.5)

    # inference arguments
    rmfdca_args = parser.add_argument_group('rmfdca_args')
    rmfdca_args.add_argument('--shrinkage_coeff', help="rmfDCA shrinkage coeff (default: 0.6)", type=float, default=0.6)
    rmfdca_args.add_argument('--uniform_pc_rate', help="uniform pseudo-count rate (default: 0.5)", type=float, default=0.5)
    rmfdca_args.add_argument('--reg_lambda_w', help="regularization coefficient for coupling parameter inference (default: 1e-4)", type=float, default=1e-4)

    args = vars(parser.parse_args(args))

    build_mrf_to_folder(**args)
