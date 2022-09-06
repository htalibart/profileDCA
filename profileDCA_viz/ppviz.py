import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from profileDCA_utils import global_variables, ppfunctions, manage_positions
from profileDCA_utils import io_management as iom
from profileDCA_align.compute_scores import *

def end_visual(tight_layout=True, show_figure=True, **kwargs):
    """ factorizing matplotlib options always used """
    if tight_layout:
        plt.tight_layout()
    plt.draw()
    if show_figure:
        plt.show()

def visualize_v_parameters(v, alphabet=global_variables.ALPHABET, start_at_1=True, tick_space=3, figsize=(10,2), center=0, **kwargs):
    xticklabels = [str(i+start_at_1) if (i%tick_space==0) else " " for i in range(0,v.shape[0])]
    v = ppfunctions.get_reordered_v(v, alphabet)
    plt.figure(figsize=figsize)
    sns.heatmap(np.transpose(v), yticklabels=alphabet, xticklabels=xticklabels, cmap="RdBu", center=center, cbar_kws={'label': r'$v_i(a)$'})
    plt.tick_params(labelsize='xx-small')
    end_visual(**kwargs)


def visualize_v_norms(v_norm, start_at_1=True, tick_space=3, figsize=(10,2), colorbar_label=r'$||v_i||$', **kwargs):
    xticklabels = [str(i+start_at_1) if (i%tick_space==0) else " " for i in range(0,len(v_norm))]
    plt.figure(figsize=figsize)
    sns.heatmap([v_norm], xticklabels=xticklabels, yticklabels=[], cmap="RdBu", center=0, cbar_kws={'label': colorbar_label})
    plt.tick_params(labelsize='xx-small')
    end_visual(**kwargs)
  

def visualize_w_norms(w_norm, start_at_1=True, tick_space=3, figsize=(10,9), tight_layout=True, colorbar_label = r'$||w_{ij}||$', **kwargs):
    xticklabels = [str(i+start_at_1) if (i%tick_space==0) else " " for i in range(0,len(w_norm))]
    plt.figure(figsize=figsize)
    sns.heatmap(w_norm, xticklabels=xticklabels, yticklabels=xticklabels, cmap="RdBu", center=0, cbar_kws={'label': colorbar_label})
    plt.tick_params(labelsize='xx-small')
    end_visual(**kwargs)



def visualize_parameters(v, v_norm, w_norm, name, alphabet=global_variables.ALPHABET, start_at_1=True, **kwargs):
    """ displays v, ||v|| and ||w|| """
    tick_space = 3
    xticklabels = [str(i+start_at_1) if (i%tick_space==0) else " " for i in range(0,v.shape[0])]

    fig, ax = plt.subplots(nrows=3, ncols=1, sharex=False, sharey=False, gridspec_kw={'height_ratios':[1,4,1]})
    plt.text(0, 2.5, '...'+name[-50:])

    v = ppfunctions.get_reordered_v(v, alphabet)
    sns.heatmap(np.transpose(v), yticklabels=alphabet, xticklabels=xticklabels, cmap="RdBu", ax=ax[0], center=0)
    ax[0].tick_params(labelsize='xx-small')
    ax[0].set_xlabel('i')
    ax[0].set_ylabel('a')
    ax[0].collections[0].colorbar.set_label("vi(a)")
    #ax[0].text(-5,10,'v',fontsize=15)

    sns.heatmap(w_norm, xticklabels=xticklabels, yticklabels=xticklabels, cmap="RdBu", center=0, ax=ax[1])
    ax[1].tick_params(labelsize='xx-small')
    ax[1].set_xlabel('i')
    ax[1].set_ylabel('j')
    ax[1].collections[0].colorbar.set_label("||wij||")

    sns.heatmap([v_norm], xticklabels=xticklabels, yticklabels=[], cmap="RdBu", center=0, ax=ax[2])
    ax[2].tick_params(labelsize='xx-small')
    ax[2].set_xlabel('i')
    ax[2].collections[0].colorbar.set_label("||vi||")

    end_visual(**kwargs)



def visualize_mrf(mrf, name="", alphabet=global_variables.ALPHABET, start_at_1=True, **kwargs):
    """ displays MRF parameters """
    visualize_parameters(mrf['v'], ppfunctions.compute_v_norms(mrf['v']), ppfunctions.compute_w_norms(mrf['w']), name, alphabet=alphabet, start_at_1=start_at_1, **kwargs)



def plot_one_vi(vi, alphabet=global_variables.ALPHABET, show_figure=True, **kwargs):
    if alphabet in kwargs:
        idx = [ALPHABET.find(a) for a in alphabet]
        new_vi = vi[idx]
    else:
        new_vi = vi
    plt.figure()
    sns.heatmap(vi.reshape(vi.shape[0], 1), square=True, cmap="RdBu", center=0, yticklabels=alphabet, xticklabels=[], **kwargs)
    plt.margins(0,0)
    plt.tight_layout()
    if show_figure:
        plt.show()


def plot_one_wij(wij, alphabet=global_variables.ALPHABET, center=0, show_figure=True, **kwargs):
    q = wij.shape[0]
    alphabet=alphabet[:q]
    wij = ppfunctions.get_reordered_wij(wij, alphabet)
    plt.figure()
    sns.heatmap(wij, cmap="RdBu", center=center, xticklabels=alphabet, yticklabels=alphabet, **kwargs)
    plt.tick_params(labelsize='xx-small')
    if show_figure:
        plt.show()

def plot_heatmap(matrix, center=0, **kwargs):
    """ plots a heatmap with seaborn """
    plt.figure()
    sns.heatmap(matrix, cmap="RdBu", center=center, **kwargs)
    plt.tick_params(labelsize='xx-small')


def visualize_v_w_scores_alignment(aligned_mrfs, aln_res_file, show_figure=True, tick_space=3, alphabet=global_variables.ALPHABET, start_at_1=False, **kwargs):

    aligned_pos = iom.get_aligned_positions_from_ppalign_output_file(aln_res_file)
    label_list = aligned_pos

    len_aln = len(aligned_pos[0])

    fig, ax = plt.subplots(nrows=4, ncols=1, sharex=False, sharey=False, gridspec_kw={'height_ratios':[1,1,1,6]})

    # v alignment : vi(a)*vk(a)
    len_v = len(aligned_mrfs[0]['v'][0])
    letter_v_scores = np.zeros((len_aln, len_v))
    for ind_i in range(len_aln):
        for a in range(len_v):
            letter_v_scores[ind_i][a] = aligned_mrfs[0]['v'][aligned_pos[0][ind_i]][a]*aligned_mrfs[1]['v'][aligned_pos[1][ind_i]][a]
    v = ppfunctions.get_reordered_v(letter_v_scores, alphabet)
    xticklabels = []
    sns.heatmap(np.transpose(v), yticklabels=alphabet, xticklabels=xticklabels, cmap="RdBu", ax=ax[0], center=0)
    ax[0].tick_params(labelsize='xx-small')


    # v scores alignment
    v_scores = get_v_scores_for_alignment(aligned_mrfs, aligned_pos, **kwargs)
    sns.heatmap([v_scores], xticklabels=[], yticklabels=['v'], cmap="RdBu", center=0, ax=ax[1])
    #ax[2].tick_params(labelsize='xx-small')

    # w scores contributions
    w_scores_sums = [sum([get_wij_wkl_score(aligned_mrfs[0]['w'][i][j],aligned_mrfs[1]['w'][k][l]) for j,l in zip(aligned_pos[0], aligned_pos[1])]) for i,k in zip(aligned_pos[0], aligned_pos[1])]
    sns.heatmap([w_scores_sums], yticklabels=['w'], xticklabels=[], cmap="RdBu", ax=ax[2], center=0)
    ax[2].tick_params(labelsize='x-small')


    # w scores
    w_scores = get_w_scores_for_alignment(aligned_mrfs, aligned_pos)
    xticklabels = [(label_list[0][k]+start_at_1,label_list[1][k]+start_at_1) for k in range(len(aligned_pos[0]))]
    xticklabels = [xi if (i%tick_space==0) else " " for i, xi in enumerate(xticklabels)]
    yticklabels=xticklabels
    sns.heatmap(w_scores, xticklabels=xticklabels, yticklabels=yticklabels, cmap="RdBu", center=0, ax=ax[3])
    ax[3].tick_params(labelsize='x-small')

    plt.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.draw()
    if show_figure:
        plt.show()



def visualize_v_scores_alignment_with_sequence(aligned_mrfs, aln_seq_with_gaps_file, show_figure=True, tick_space=3, alphabet=global_variables.ALPHABET, start_at_1=False, **kwargs):

    aligned_seq_pos = iom.get_aligned_positions_with_gaps_from_ppalign_output_file(aln_seq_with_gaps_file)
    seq_pos_to_mrf_pos_lists = [manage_positions.reverse_pos1_to_pos2(mrf['mrf_pos_to_seq_pos'], len(mrf['v_full'])) for mrf in aligned_mrfs]
    
    label_list = aligned_seq_pos

    len_seq_aln = len(aligned_seq_pos[0])

    #fig, ax = plt.subplots(nrows=5, ncols=1, sharex=False, sharey=False, gridspec_kw={'height_ratios':[1,1,1,1,6]})
    fig, ax = plt.subplots(nrows=3,ncols=1, sharex=False, sharey=False, gridspec_kw={'height_ratios':[1,1,1]})
    
    # v full alignment : vi(a)*vk(a) when aligned, 10 when gap in A, -10 when gap in B
    len_v_full = len(aligned_mrfs[0]['v_full'][0])
    letter_v_scores = np.zeros((len_seq_aln, len_v_full))
    for ind_i in range(len_seq_aln):
        for a in range(len_v_full):
            pos_in_A = aligned_seq_pos[0][ind_i]
            if pos_in_A!='-':
                pos_in_B = aligned_seq_pos[1][ind_i]
                if pos_in_B!='-':
                    letter_v_scores[ind_i][a] = aligned_mrfs[0]['v_full'][pos_in_A][a]*aligned_mrfs[1]['v_full'][pos_in_B][a]
                else:
                    letter_v_scores[ind_i][a] = -1
            else:
                letter_v_scores[ind_i][a] = 1
    v_full = ppfunctions.get_reordered_v(letter_v_scores, alphabet)
    xticklabels = []
    sns.heatmap(np.transpose(v_full), yticklabels=alphabet, xticklabels=xticklabels, cmap="BrBG", ax=ax[0], center=0)
    ax[0].tick_params(labelsize='xx-small')
    
    # v alignment : vi(a)*vk(b) when aligned and in MRF, 0 otherwise
    len_v = len(aligned_mrfs[0]['v'][0])
    letter_v_scores = np.zeros((len_seq_aln, len_v))
    v_scores = np.zeros((len_seq_aln))
    for ind_i in range(len_seq_aln):
        pos_in_A = aligned_seq_pos[0][ind_i]
        pos_in_B = aligned_seq_pos[1][ind_i]
        if (pos_in_A!='-') and (pos_in_B!='-'):
            pos_in_mrf_A = seq_pos_to_mrf_pos_lists[0][pos_in_A]
            pos_in_mrf_B = seq_pos_to_mrf_pos_lists[1][pos_in_B]
            if (pos_in_mrf_A!=None) and (pos_in_mrf_B!=None):
                for a in range(len_v):
                    letter_v_scores[ind_i][a] = aligned_mrfs[0]['v'][pos_in_mrf_A][a]*aligned_mrfs[1]['v'][pos_in_mrf_B][a]
                v_scores[ind_i] = get_vi_vk_score(aligned_mrfs[0]['v'][pos_in_mrf_A], aligned_mrfs[1]['v'][pos_in_mrf_B])
    v = ppfunctions.get_reordered_v(letter_v_scores, alphabet)
    xticklabels = []
    sns.heatmap(np.transpose(v), yticklabels=alphabet, xticklabels=xticklabels, cmap="RdBu", ax=ax[1], center=0)
    ax[1].tick_params(labelsize='xx-small')

    # v scores alignment
    sns.heatmap([v_scores], xticklabels=[], yticklabels=['v'], cmap="RdBu", center=0, ax=ax[2])

    plt.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.draw()
    if show_figure:
        plt.show()
