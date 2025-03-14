import numpy as np
from profileDCA_build import msa_statistics, trim, pseudocounts, covariance_processing
from profileDCA_utils import io_management as iom
from profileDCA_utils import global_variables

def single_frequencies_to_fields(single_freqs):
    """ inputs single frequency array @single_freqs of shape Lx20,
        outputs field array directly blindly computed from frequencies
    """
    L = single_freqs.shape[0]
    assert(single_freqs.shape[1]==20)
    lf = np.log(single_freqs)
    v_star = lf - np.tile(np.mean(lf, axis=1), (single_freqs.shape[1],1)).T
    return v_star

def compute_fields(fi_with_gaps, uniform_pc_rate=0.5):
    """ input Lx21 array of single frequencies with gaps @fi_with_gaps and a desired uniform pseudo-count rate @uniform_pc_rate 
        outputs field array for profileDCA, i.e. a Lx20 array of fields where gap symbols are ignored and uniform pseudo-counts are applied
    """
    fi = msa_statistics.fi_with_gaps_to_fi_without_gaps(fi_with_gaps)
    fi_pc = pseudocounts.apply_uniform_pseudocounts_to_single_frequencies(fi, uniform_pc_rate)
    v = single_frequencies_to_fields(fi_pc)
    return v

def J_to_w(J):
    """ transforms matrix (L*19,L*19) matrix J to a (L,L,20,20) coupling array by unravelling it and filling values for the 20^th amino acid """
    L = J.shape[0]//19
    w = np.zeros((L,L,20,20))
    w[:,:,:19,:19] = np.transpose(J.reshape((L,19,L,19)),axes=(0,2,1,3))
    for i in range(L-1):
        for j in range(i+1,L):
            w[i,j,19,:19] = -np.sum(w[i,j,:19,:19], axis=0)
            w[i,j,:,19] = -np.sum(w[i,j,:,:19], axis=1)
            w[j,i] = np.transpose(w[i,j])
    w[range(L),range(L),:,:]=0
    return w


def regularized_inversion(C, reg_lambda_w):
    """ returns the solution of min_J ||JC+I||_2 + lw||J||_2 """
    C2 = np.matmul(C,C)
    I = np.identity(C.shape[0])
    Cinv = np.linalg.inv(C2+reg_lambda_w*I)
    return -np.matmul(C, Cinv)

def compute_couplings(fi, fij, covariance_matrix_norm_threshold, uniform_pc_rate, shrinkage_coeff, reg_lambda_w):
    """ inputs arrays of frequencies computed with gaps @fi (L,21) and @fij (L,L,21,21),
        minimum threshold for each covariance matrix block ||Cij||>@covariance_matrix_norm_threshold,
        desired uniform pseudo-count rate @uniform_pc_rate,
        desired uniform covariance matrix shrinkage coefficient @shrinkage_coeff,
        desired regularization coefficient @reg_lambda_w
        outputs coupling array
    """
    C = covariance_processing.compute_covariance_matrix_without_gaps(fi, fij)
    C = covariance_processing.apply_norm_threshold(C, covariance_matrix_norm_threshold)
    fi = msa_statistics.fi_with_gaps_to_fi_without_gaps(fi)
    fi_pc = pseudocounts.apply_uniform_pseudocounts_to_single_frequencies(fi, uniform_pc_rate)
    C = pseudocounts.apply_uniform_pseudocounts_to_covariance_matrix(C, fi_pc, uniform_pc_rate)
    C = covariance_processing.shrink_covariance_matrix_towards_uniform(C, shrinkage_coeff)
    Czero = covariance_processing.prepare_covariance_matrix_for_zero_sum_gauge(C)
    J = regularized_inversion(Czero, reg_lambda_w)
    w = J_to_w(J)
    return w


def msa_array_to_mrf(msa_int_array, max_gap_v=1, max_gap_w=0.9, uniform_pc_rate=0.5, apply_covariance_matrix_threshold=True, covariance_matrix_norm_threshold=None, shrinkage_coeff=0.7, reg_lambda_w=1e-4, alphabet_dict=global_variables.ALPHABET_DICT):
    """
        inputs MSA as an integer array @msa_int_array and outputs a MRF representing this MSA in dictionary format {'v': array of fields, 'w': array of couplings, etc.}
        @max_gap_v (resp. @max_gap_w) is the maximum fraction of gaps allowed in a column so that it is included in fields (resp. couplings)
        @uniform_pc_rate is the uniform pseudo-count rate
        when @apply_covariance_matrix_threshold is True, covariance matrix is processed so that blocks (i,j) where ||Cij||<@covariance_matrix_norm_threshold are set to 0
        when @covariance_matrix_norm_threshold is None, this threshold is computed as a function of N
        @shrinkage_coeff is the coefficient for covariance matrix shrinkage towards uniform distribution
        @reg_lambda_w is the coefficient for L2 regularization
    """
    mrf = {}

    L = msa_int_array.shape[1]

    fi_whole_msa = msa_statistics.compute_single_frequencies(msa_int_array, alphabet_dict=alphabet_dict)
    mrf['v_full'] = compute_fields(fi_whole_msa, uniform_pc_rate=uniform_pc_rate)
   
    msa_trimmed_for_v, mrf['mrf_pos_to_seq_pos'] = trim.trim_int_msa_array(msa_int_array, max_gap_v, alphabet_dict['-'])
    mrf['v'] = mrf['v_full'][mrf['mrf_pos_to_seq_pos'],:]
    
    L_v = msa_trimmed_for_v.shape[1]
    mrf['w'] = np.zeros((L_v,L_v,20,20))
    msa_trimmed_for_w, pos_w_trim = trim.trim_int_msa_array(msa_trimmed_for_v, max_gap_w, alphabet_dict['-'])
    N, L_trimmed = msa_trimmed_for_w.shape
    if L_trimmed>1:
        fi_trimmed = fi_whole_msa[mrf['mrf_pos_to_seq_pos'],:][pos_w_trim,:]
        fij_trimmed = msa_statistics.compute_double_frequencies(msa_trimmed_for_w, alphabet_dict=alphabet_dict)
        if apply_covariance_matrix_threshold:
            if covariance_matrix_norm_threshold==None:
                covariance_matrix_norm_threshold = covariance_processing.get_norm_covariance_matrix_threshold(N)
        else:
            covariance_matrix_norm_threshold = 0
        w_trim = compute_couplings(fi_trimmed, fij_trimmed, covariance_matrix_norm_threshold, uniform_pc_rate, shrinkage_coeff, reg_lambda_w)
        mrf['w'][np.ix_(pos_w_trim, pos_w_trim)] = w_trim
    return mrf
