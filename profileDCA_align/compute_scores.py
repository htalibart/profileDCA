import numpy as np
from profileDCA_build import pseudocounts, mrf_inference
from profileDCA_utils import global_variables
from profileDCA_utils import io_management as iom

def get_vi_vk_score(vi, vk, offset_v=0):
    """ returns the similarity score of field vectors @vi and @vk """
    return np.dot(vi,vk)-offset_v


def compute_v_scores(mrf1, mrf2, v_score_function, offset_v, **kwargs):
    """ returns an np.array @v_scores where v_scores[i][k] is the similarity score between field vi in Potts model @mrf1 and field vk in Potts model @mrf2 """
    v_scores = np.zeros((mrf1['v'].shape[0], mrf2['v'].shape[0]))
    for i in range(mrf1['v'].shape[0]):
        for k in range(mrf2['v'].shape[0]):
            v_scores[i][k] = get_vi_vk_score(mrf1['v'][i], mrf2['v'][k], offset_v=offset_v)
    return v_scores


def get_wij_wkl_score(wij, wkl):
    """ returns the similarity score of couplint matrices @wij and @wkl """
    return np.dot(wij.flatten(), wkl.flatten())


def get_epsilon(epsilon_sim, selfscores):
    """ gives epsilon for s(A,B) so that 2*s(A,B)/(s(A,A)+s(B,B)) < @epsilon_sim (where s(A,A) and s(B,B) are @selfscores)"""
    return (epsilon_sim/2)*sum(selfscores)


def compute_self_w_scores(mrf, edges_map):
    """ @w_score[i,j] is score for aligning w[i,j] with w[i,j] in Potts model @mrf for each considered edge in @edges_map """
    w_score = 0
    L = mrf['v'].shape[0]
    for i in range(L-1):
        for j in range(i+1,L):
            if edges_map[i][j]:
                w_score+=get_wij_wkl_score(mrf['w'][i][j], mrf['w'][i][j])
    return w_score


def compute_selfscore(mrf, edges_map, alpha_w=1, offset_v=0, use_v=True, use_w=True):
    """ score for aligning @mrf with itself """
    if use_v:
        v_score = sum([get_vi_vk_score(vi, vi, offset_v=offset_v) for vi in mrf['v']])
    else:
        v_score = 0
    if use_w:
        w_score = compute_self_w_scores(mrf, edges_map)
    else:
        w_score = 0
    selfcomp = v_score+alpha_w*w_score
    return selfcomp

def get_v_scores_for_alignment(aligned_potts_models, aligned_pos, offset_v=0, **kwargs):
    """ @v_scores[i] is the score for the alignment of field vectors at positions in the two @aligned_potts_models that are aligned at position i in the PPalign alignment given by @aligned_positions_dict """
    v_scores = np.array([get_vi_vk_score(aligned_potts_models[0]['v'][i],
                                         aligned_potts_models[1]['v'][k],
                                         offset_v=offset_v) for i,k in zip(aligned_pos[0], aligned_pos[1])])
    return v_scores


def get_w_scores_for_alignment_given_w(w_list, aligned_pos):
    """ @w_scores[i,j] is the score for the alignment of couplings from the two @aligned_potts_models that are aligned at position i and j in the PPalign alignment given by @aligned_positions_dict """
    L = len(aligned_pos[0])
    w_scores = np.zeros((L,L))
    for ind_i in range(L-1):
        for ind_j in range(ind_i+1,L):
            w_scores[ind_i,ind_j] = get_wij_wkl_score(w_list[0][aligned_pos[0][ind_i],aligned_pos[0][ind_j]], w_list[1][aligned_pos[1][ind_i],aligned_pos[1][ind_j]])
            w_scores[ind_j,ind_i]=w_scores[ind_i,ind_j]
    return w_scores


def get_w_scores_for_alignment(aligned_potts_models, aligned_pos):
    return get_w_scores_for_alignment_given_w([mrf['w'] for mrf in aligned_potts_models], aligned_pos)


def get_offset_for_uniform_pc_rate(uniform_pc_rate):
    bg_fi = np.array(global_variables.AA_BACKGROUND_FREQUENCIES).reshape((1,20))
    bg_fi_pc = pseudocounts.apply_uniform_pseudocounts_to_single_frequencies(bg_fi, uniform_pc_rate)
    bg_vi_pc = mrf_inference.single_frequencies_to_fields(bg_fi_pc)[0]
    return get_vi_vk_score(bg_vi_pc, bg_vi_pc, offset_v=0) 


def get_offset(potts_folders):
    """ computes offset from README files in Potts folders. Raises Exception if no README or if uniform_pc_rate is not the same """
    uniform_pcs = []
    for pf in potts_folders:
        readme_file = pf/"README.txt"
        if not readme_file.is_file():
            raise Exception("No README.txt in "+str(pf)+". Offset should be set manually with offset_v option")
        inference_params = iom.get_parameters_from_readme_file(readme_file)
        uniform_pcs.append(float(inference_params['uniform_pc_rate']))
    if not uniform_pcs[0]==uniform_pcs[1]:
        raise Exception("Potts models were trained with a different uniform pseudo-count rate... If you still want to align them, offset_v should be set manually")
    else:
        uniform_pc_rate = uniform_pcs[0]
    return get_offset_for_uniform_pc_rate(uniform_pc_rate)
