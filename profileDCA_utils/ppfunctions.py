import numpy as np
from Bio import pairwise2

from profileDCA_utils import global_variables

def get_letter_index(letter, alphabet_dict):
    """ inputs dictionary @alphabet_dict {letter: index}
        outputs index if letter is in dictionary, raises exception otherwise """
    if letter in alphabet_dict:
        return alphabet_dict[letter]
    else:
        raise Exception("letter "+str(letter)+" not in alphabet")

def compute_w_norms(w):
    """ outputs (LxL) array of Frobenius norms for coupling parameters @w """
    L = w.shape[0]
    w_norms = np.zeros((L,L))
    for i in range(0, L-1):
        for j in range(i+1, L):
            w_norms[i][j] = np.linalg.norm(w[i,j])
            w_norms[j][i] = w_norms[i][j]
    return w_norms

def compute_v_norms(v):
    """ outputs list of L2 norms for field parameters @v """
    return np.array([np.linalg.norm(vi) for vi in v])


def get_reordered_v(v, alphabet_to, alphabet_from=global_variables.ALPHABET):
    """ reorders all vi from alphabet @alphabet_to to alphabet @alphabet_from """
    q = v.shape[1]
    idx = [alphabet_from[:q].find(a) for a in alphabet_to[:q]]
    return v[:,idx]


def get_reordered_wij(wij, alphabet_to, alphabet_from=global_variables.ALPHABET):
    """ reorders all wij from alphabet @alphabet_to to alphabet @alphabet_from """
    q = wij.shape[0]
    idx = [alphabet_from[:q].find(a) for a in alphabet_to[:q]]
    new_wij = np.zeros_like(wij)
    for a in range(len(alphabet_to)):
        for b in range(len(alphabet_to)):
            new_wij[a][b] = wij[idx[a]][idx[b]]
    return new_wij

def get_reordered_w(w, alphabet_to, alphabet_from=global_variables.ALPHABET):
    """ reorders all wij of @w from alphabet @alphabet_to to alphabet @alphabet_from """
    new_w = np.zeros_like(w)
    L = w.shape[0]
    for i in range(L-1):
        for j in range(i+1,L):
            new_w[i,j] = get_reordered_wij(w[i,j], alphabet_to=alphabet_to, alphabet_from=alphabet_from)
            new_w[j,i] = np.transpose(new_w[i,j])
    return new_w

def get_pos_first_seq_to_second_seq(first_seq, second_seq):
    """ returns dictionary d[pos_in_first_seq] = pos_in_second_seq """
    gap_char='.'
    alns = pairwise2.align.globalxx(first_seq, second_seq, gap_char=gap_char)
    top_aln = alns[0]
    aln_first, aln_second, score, begin, end = top_aln
    first_pos = 0
    second_pos = 0
    pos_dict_first_seq_to_second_seq = [None for k in range(len(first_seq))]
    for i in range(len(aln_first)):
        if (aln_first[i]==gap_char) and (aln_second[i]!=gap_char):
            second_pos+=1
        elif (aln_first[i]!=gap_char) and (aln_second[i]==gap_char):
            pos_dict_first_seq_to_second_seq[first_pos] = None
            first_pos+=1
        elif (aln_first[i]!=gap_char) and (aln_second[i]!=gap_char):
            pos_dict_first_seq_to_second_seq[first_pos] = second_pos
            first_pos+=1
            second_pos+=1
    return pos_dict_first_seq_to_second_seq

def unravel_covariance(C):
    """ inputs (L*(q-1),L*(q-1)) array, outputs Lx(q-1)xLx(q-1) array, where q=20 """
    L = C.shape[0]//19
    return np.transpose(C.reshape((L,19,L,19)),axes=(0,2,1,3))
