import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from profileDCA_utils import global_variables, ppfunctions

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

