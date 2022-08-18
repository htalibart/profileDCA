import ctypes, time, math, pathlib, tempfile, shutil
import pandas as pd
import numpy as np

from profileDCA_utils import ppfunctions, global_variables, manage_positions
from profileDCA_utils import io_management as iom
from profileDCA_align.compute_scores import *
from profileDCA_build import __main__ as ppbuild_main

import pkg_resources
PPALIGN_CPP_LIBRARY = pkg_resources.resource_filename('profileDCA_align', 'ppalign_solver.so')
PPALIGN_SOLVER = ctypes.CDLL(PPALIGN_CPP_LIBRARY)

INFINITY = 1000000000

def get_edges_map(w): #TODO w percent
    edges_map = 1*(ppfunctions.compute_w_norms(w)>0)
    for i in range(w.shape[0]):
        edges_map[i,i] = 0
    return edges_map

def align_two_mrfs_with_solver(mrfs, output_folder, n_limit_param=INFINITY, iter_limit_param=1000, t_limit=36000, disp_level=1, epsilon_sim=0.001, use_w=True, use_v=True, gamma=1.0, theta=0.9, stepsize_min=0.000000005, nb_non_increasing_steps_max=500, alpha_w=1, sim_min=-100, offset_v=0, gap_open=8, gap_extend=0, insert_opens=None, insert_extends=None, costly_end_gap_left=False, costly_end_gap_right=False, remove_negative_couplings=False, **kwargs):

    mrf_lengths = [mrf['v'].shape[0] for mrf in mrfs]
    if not all([mrf_length>0 for mrf_length in mrf_lengths]):
        if not any(mrf_lengths):
            aligned_positions_with_gaps = [[],[]]
        elif (mrf_lengths[0]==0):
            aligned_positions_with_gaps = [['-' for k in range(mrf_lengths[1])], [k for k in range(mrf_lengths[1])]]
        else:
            aligned_positions_with_gaps = [[k for k in range(mrf_lengths[0])], ['-' for k in range(mrf_lengths[0])]]
        return {
                'aligned_positions':[[],[]],
                'aligned_positions_with_gaps': aligned_positions_with_gaps,
                'infos_solver': None
                }

    else:
        # handle output files and folder
        if not output_folder.is_dir():
            iom.create_folder(output_folder)
        aln_res_file = output_folder/"aln.csv"
        info_res_file = output_folder/"info.csv"
        aln_with_gaps_res_file = output_folder/"aln_with_gaps.csv"


       # DEFINE TYPES FOR CTYPES
        c_float_p = ctypes.POINTER(ctypes.c_float) # pointer to float
        c_int_p = ctypes.POINTER(ctypes.c_int) # pointer to int


        # FIELDS
        if use_v:
            vs = [mrf['v'] for mrf in mrfs]
            v_flats = [np.ascontiguousarray(v.flatten()) for v in vs] # flatten for ctypes
        else:
            v_flats = [np.ascontiguousarray(np.zeros((mrf['v'].shape[0],20))).flatten() for mrf in mrfs]
        c_vs = [vflat.astype(np.float32).ctypes.data_as(c_float_p) for vflat in v_flats]


        # COUPLINGS
        if use_w:
            edges_maps = [get_edges_map(mrf['w']) for mrf in mrfs]
            # always give C++ program coupling matrices 20x20 => add zeroes if necessary
            w_flats = [np.ascontiguousarray(mrf['w'].flatten()) for mrf in mrfs]
        else: # if no coupling: edge map where edges are all 0
            edges_maps = [np.zeros((mrf['w'].shape[0:2])) for mrf in mrfs]
            w_flats = [np.ascontiguousarray(np.zeros((mrf['w'].shape[0],mrf['w'].shape[1],20,20))).flatten() for mrf in mrfs]
        c_ws = [wflat.astype(np.float32).ctypes.data_as(c_float_p) for wflat in w_flats]

        selfcomps = [compute_selfscore(mrf, edges_map, alpha_w=alpha_w, offset_v=offset_v, use_v=use_v, use_w=use_w) for mrf, edges_map in zip(mrfs, edges_maps)] #s(A,A) and s(B,B)
        epsilon = get_epsilon(epsilon_sim, selfcomps) # epsilon so that 2*s(A,B)/(s(A,A)+s(B,B)) < @epsilon_sim
        score_min = (1/2)*sim_min*sum(selfcomps); # stop computation if LB < score_min 

        c_edges_maps = [np.ascontiguousarray(edges_map.flatten(), dtype=np.int32).ctypes.data_as(c_int_p) for edges_map in edges_maps]

        # INSERTIONS
        if insert_opens==None:
            insert_opens = [np.ascontiguousarray(np.hstack((np.zeros((1)), np.ones((mrf_lengths[mrf_ind]-1))*gap_open, np.zeros(1)))) for mrf_ind in range(2)]
        if insert_extends==None:
            insert_extends = [np.ascontiguousarray(np.hstack((np.zeros((1)), np.ones((mrf_lengths[mrf_ind]-1))*gap_extend, np.zeros(1)))) for mrf_ind in range(2)]
        if costly_end_gap_left:
            for mrf_ind in range(2):
                insert_opens[mrf_ind][0] = gap_open
                insert_extends[mrf_ind][0] = gap_extend
        if costly_end_gap_right:
            for mrf_ind in range(2):
                insert_opens[mrf_ind][-1] = gap_open
                insert_extends[mrf_ind][-1] = gap_extend
        c_insert_opens = [insert_open.astype(np.float32).ctypes.data_as(c_float_p) for insert_open in insert_opens]
        c_insert_extends = [insert_extend.astype(np.float32).ctypes.data_as(c_float_p) for insert_extend in insert_extends]

        PPALIGN_SOLVER.call_from_python.argtypes=[c_float_p, # c_vs[0]
                c_float_p, # c_vs[1]
                c_float_p, # c_ws[0]
                c_float_p, # c_ws[1]
                ctypes.c_int, # LA
                ctypes.c_int, # LB
                c_int_p, # edge map A
                c_int_p, # edge map B
                ctypes.c_double, # self alignemnt score A
                ctypes.c_double, # self alignment score B
                c_float_p, # insertion penalties open A
                c_float_p, # insertion penalties open B
                c_float_p, # insertion penalties extend A
                c_float_p, # insertion penalties extend B
                ctypes.c_char_p, # aln_res_file
                ctypes.c_char_p, # info_res_file
                ctypes.c_char_p, # aln_with_gaps_res_file
                ctypes.c_int, # n_limit_param
                ctypes.c_int, # iter_limit_param
                ctypes.c_double, # t_limit
                ctypes.c_int, # disp level
                ctypes.c_double, # epsilon
                ctypes.c_double, # gamma
                ctypes.c_double, # theta
                ctypes.c_double, # stepsize_min
                ctypes.c_int, # nb_non_increasing_steps_max
                ctypes.c_double,
                ctypes.c_double,
                ctypes.c_double]


        PPALIGN_SOLVER.call_from_python(
                *c_vs,
                *c_ws,
                *[ctypes.c_int(mrf_length) for mrf_length in mrf_lengths],
                *c_edges_maps,
                *[ctypes.c_double(selfcomp) for selfcomp in selfcomps],
                *c_insert_opens,
                *c_insert_extends,
                ctypes.c_char_p(str(aln_res_file).encode('utf-8')),
                ctypes.c_char_p(str(info_res_file).encode('utf-8')),
                ctypes.c_char_p(str(aln_with_gaps_res_file).encode('utf-8')),
                ctypes.c_int(n_limit_param),
                ctypes.c_int(iter_limit_param),
                ctypes.c_double(t_limit),
                ctypes.c_int(disp_level),
                ctypes.c_double(epsilon),
                ctypes.c_double(gamma),
                ctypes.c_double(theta),
                ctypes.c_double(stepsize_min),
                ctypes.c_int(nb_non_increasing_steps_max),
                ctypes.c_double(score_min),
                ctypes.c_double(alpha_w),
                ctypes.c_double(offset_v)
                )

        infos_solver = iom.get_infos_solver_dict_from_ppalign_output_file(info_res_file)

        if not math.isnan(infos_solver['similarity_global']):
            aligned_positions = iom.get_aligned_positions_from_ppalign_output_file(aln_res_file)
            aligned_positions_with_gaps = iom.get_aligned_positions_with_gaps_from_ppalign_output_file(aln_with_gaps_res_file)
        else:
            aligned_positions = []
            aligned_positions_with_gaps = []

        return {
                'aligned_positions':aligned_positions,
                'aligned_positions_with_gaps': aligned_positions_with_gaps,
                'infos_solver': infos_solver
                }




def align_missing_positions(aligned_sequence_positions, mrfs, pc_tau=0.5, alphabet=global_variables.ALPHABET, **kwargs):
    kwargs['use_w']=False
    previously_aligned = [-1,-1]
    whole_alignment = [[],[]]
    for current_aln_pos in range(len(aligned_sequence_positions[0])):
        just_aligned = [aligned_sequence_positions[seq_index][current_aln_pos] for seq_index in range(2)]
        to_be_aligned = [[previously_aligned[seq_index]+1,just_aligned[seq_index]-1] for seq_index in range(2)]
        sub_profiles = [ppbuild_main.get_sub_v(mrfs[seq_index], *to_be_aligned[seq_index]) for seq_index in range(2)]
        tmp_output_folder = pathlib.Path(tempfile.mkdtemp())
        if previously_aligned==[-1,-1]:
            costly_end_gap_left=False
        else:
            costly_end_gap_left=True
        res_subalignment = align_two_mrfs_with_solver(sub_profiles, tmp_output_folder, disp_level=1, costly_end_gap_left=costly_end_gap_left, costly_end_gap_right=True, **kwargs)
        shutil.rmtree(tmp_output_folder)
        for seq_index in range(2):
            sub_aligned_pos = res_subalignment['aligned_positions_with_gaps'][seq_index]
            for sub_aln_pos in sub_aligned_pos:
                if sub_aln_pos!='-':
                    pos_to_add = sub_aln_pos+to_be_aligned[seq_index][0]
                else:
                    pos_to_add='-'
                whole_alignment[seq_index].append(pos_to_add)
            whole_alignment[seq_index].append(just_aligned[seq_index])
        previously_aligned = just_aligned
    # align sequence ends
    to_be_aligned = [[previously_aligned[seq_index]+1,len(mrfs[seq_index]['v_full'])-1] for seq_index in range(2)]
    sub_profiles = [ppbuild_main.get_sub_v(mrfs[seq_index], *to_be_aligned[seq_index]) for seq_index in range(2)]
    tmp_output_folder = pathlib.Path(tempfile.mkdtemp())
    res_subalignment = align_two_mrfs_with_solver(sub_profiles, tmp_output_folder, disp_level=1, costly_end_gap_left=True, costly_end_gap_right=False, **kwargs)
    shutil.rmtree(tmp_output_folder)
    for seq_index in range(2):
        sub_aligned_pos = res_subalignment['aligned_positions_with_gaps'][seq_index]
        for sub_aln_pos in sub_aligned_pos:
            if sub_aln_pos!='-':
                pos_to_add = sub_aln_pos+to_be_aligned[seq_index][0]
            else:
                pos_to_add='-'
            whole_alignment[seq_index].append(pos_to_add)
    return whole_alignment




def align_two_mrfs_full(mrfs, output_folder, gap_open=8, **kwargs):
    # total PPalign time starts now
    time_start = time.time()
      
    insert_opens = [np.array(manage_positions.get_insert_opens_for_mrf_pos_to_seq_pos(gap_open, mrf['mrf_pos_to_seq_pos'])) for mrf in mrfs]
    res_dict = align_two_mrfs_with_solver(mrfs, output_folder, insert_opens=insert_opens, **kwargs)
    # if Potts models were successfully aligned
    if len(res_dict['aligned_positions'])>0:
        res_dict['aligned_sequence_positions_confident'] = manage_positions.translate_aligned_positions_without_gaps(res_dict['aligned_positions'], [mrf['mrf_pos_to_seq_pos'] for mrf in mrfs])
        res_dict['aligned_sequence_positions_with_gaps'] = align_missing_positions(res_dict['aligned_sequence_positions_confident'], mrfs, gap_open=gap_open, **kwargs)
    else:
        res_dict['aligned_sequence_positions_with_gaps'] = []

    total_computation_time = time.time()-time_start
    res_dict['infos_solver']['total_time']=total_computation_time

    return res_dict

