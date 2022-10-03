import argparse, sys, time, pathlib

from profileDCA_utils import io_management as iom
from profileDCA_utils import manage_positions
from profileDCA_align.potts_model_alignment import *

def align_objects_and_handle_files(mrfs, output_folder, sequence_files, **kwargs):
    # write readme with all input arguments used
    iom.write_readme(output_folder, **kwargs)

    # get aligned sequences
    res_dict = align_two_mrfs_full(mrfs, output_folder, **kwargs)
    print("Total time : "+str(res_dict['infos_solver']["total_time"]))

    # if Potts models were successfully aligned
    if len(res_dict['aligned_sequence_positions_with_gaps'])>0:
        # write aligned sequence positions in csv file
        iom.write_positions_to_csv(res_dict['aligned_sequence_positions_with_gaps'], output_folder/("aln_sequences_with_gaps.csv"))

        # write aligned sequences in fasta file
        output_fasta_file = output_folder/("aligned_sequences.fasta")
        sequences = []
        names = []
        for seq_file in sequence_files:
            seq, name = iom.get_first_record_sequence_and_name(seq_file)
            sequences.append(seq)
            names.append(name)
        manage_positions.write_aligned_sequences_in_fasta_file_with_lowercase(res_dict['aligned_sequence_positions_with_gaps'], res_dict['aligned_sequence_positions_confident'], sequences, names, output_fasta_file)

    info_res_file = output_folder/"info.csv"
    iom.dict_to_csv(res_dict['infos_solver'], info_res_file)
    return res_dict


def main(args=sys.argv[1:]):
    parser = argparse.ArgumentParser()
    parser.add_argument('-f1', '--potts_folder_1', help="Folder containing files for sequence 1", type=pathlib.Path)
    parser.add_argument('-f2', '--potts_folder_2', help="Folder containing files for sequence 2", type=pathlib.Path)
    parser.add_argument('-o', '--output_folder', help="Output folder (if not specified : output_ppalign/[DATE]/)", type=pathlib.Path)


    # options alignement
    parser.add_argument('-nw', '--no_w', help="Don't use w scores (default : False)", action='store_true')
    parser.add_argument('-nv', '--no_v', help="Don't use v scores (default : False)", action='store_true')
    parser.add_argument('--alpha_w', help="coefficient before w score", default=1, type=float)
    parser.add_argument('--offset_v', help="score offset for v parameters", default=1, type=float)
    parser.add_argument('-go', '--gap_open', help="gap open", default=14, type=float)
    parser.add_argument('-ge', '--gap_extend', help="gap extend", default=0, type=float)
    parser.add_argument('--no_free_end_gaps', help="End gaps are not free (default: they are)", action='store_true')

    # solver options
    parser.add_argument('-t', '--t_limit', help="solver : time limit in seconds (default : 36000)", type=float, default=36000)
    parser.add_argument('--epsilon_sim', help="solver : max 2*(UB-LB)/(s(A,A)+s(B,B)) (default : 0.001)", default=0.001, type=float)
    parser.add_argument('-sim_min', '--sim_min', help='if similarity score is below sim_min, solver considers that the Potts models are not similar and stops. (default : -100)', type=float, default=-100)

    args = vars(parser.parse_args(args))

    if args["gap_extend"]>args["gap_open"]:
        raise Exception("Gap extend must be smaller than gap open (convergence issues with the solver)")


    # CREATE_FOLDER IF NOT EXISTING
    if args["output_folder"] is None:
        general_output_folder = iom.create_folder("output_ppalign")
        output_folder = general_output_folder / time.strftime("%Y%m%d-%H%M%S")
    else:
        output_folder = args["output_folder"]
    iom.create_folder(output_folder)
    del args["output_folder"]


    # HANDLING NO W / NO V ARGUMENTS
    if args["no_w"]:
        args["w_threshold"]=float('inf')
        args["use_w"]=False
    else:
        args["use_w"]=True
    del args["no_w"]

    args["use_v"]= not args["no_v"]
    del args["no_v"]



    # GET MRFS FROM FOLDERS
    mrfs = []
    sequence_files = []
    for k in range(1,3):
        pf = args["potts_folder_"+str(k)]
        if pf is not None:
            mrfs.append(iom.mrf_from_folder(pf))
            sequence_files.append(pf/"sequence.fasta")
    
    align_objects_and_handle_files(mrfs, output_folder, sequence_files, **args)



if __name__=="__main__":
    main()
