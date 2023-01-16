import numpy as np
from profileDCA_utils.global_variables import ALPHABET, ALPHABET_DICT


def compute_single_frequencies(msa, alphabet_dict=ALPHABET_DICT):
    """ input: MSA (np int array object) @msa, dictionary @alphabet_dict where alphabet_dict[a] is the index of letter a
        output: single frequencies in (L,q) np array """
    N,L = msa.shape
    q = len(alphabet_dict)
    single_counts = np.zeros((L,q))
    for i in range(L):
        unique, counts = np.unique(msa[:,i], return_counts=True)
        unique = unique.astype(int)
        single_counts[i,unique]=counts
    single_f = single_counts/N
    return single_f
                

def compute_double_frequencies(msa, alphabet_dict=ALPHABET_DICT):
    """ input: MSA (np int array object) @msa, dictionary @alphabet_dict where alphabet_dict[a] is the index of letter a
        output: double frequencies in (L,L,q,q) np array """
    N,L = msa.shape
    q = len(alphabet_dict)
    double_counts = np.zeros((L,L,q,q))
    for i in range(L):
        for j in range(i,L):
            for n in range(N):
                if (msa[n,i] in alphabet_dict.values()) and (msa[n,j] in alphabet_dict.values()):
                    double_counts[i,j,msa[n,i],msa[n,j]]+=1
            double_counts[j,i] = np.transpose(double_counts[i,j])
    double_f = double_counts/N
    return double_f


def fi_with_gaps_to_fi_without_gaps(fi):
    """ inputs frequency array @fi of shape Lx20
        outputs frequency array of shape Lx19 where frequencies are computed without gap symbol
    """
    fiq = np.tile(fi[:,20],(20,1)).T
    return 1/(1-fiq)*fi[:,:20]

def compute_frequencies(msa, alphabet_dict=ALPHABET_DICT):
    """ input: MSA (np int array object) @msa, dictionary @alphabet_dict where alphabet_dict[a] is the index of letter a
        output: single frequencies in (L,q) np array, double frequencies in (L,L,q,q) np array """
    return compute_single_frequencies(msa, alphabet_dict=alphabet_dict), compute_double_frequencies(msa, alphabet_dict=ALPHABET_DICT)
