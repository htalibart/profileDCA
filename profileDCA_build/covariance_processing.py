import numpy as np
from scipy.linalg import block_diag

def compute_gap_correction_matrix(fi, fij):
    L = fi.shape[0]
    return 1/(np.ones((L,L))-np.add.outer(fi[:,20],fi[:,20])+fij[:,:,20,20])


def compute_covariance_matrix_without_gaps(fi_with_gaps, fij_with_gaps):
    assert(fi_with_gaps.shape[1]==21)
    assert(fij_with_gaps.shape[2]==21)
    L = fi_with_gaps.shape[0]
    gap_correction_matrix = np.kron(compute_gap_correction_matrix(fi_with_gaps, fij_with_gaps), np.ones((19,19))) # Kronecker product: 20x20 -> (Lx20)x(Lx20)
    fi_j = fi_with_gaps[:,np.newaxis,:19]-fij_with_gaps[:,:,:19,20] # LxLx20 fi_j(a) = fi(a)-fij(a,-)
    fj_i = fi_with_gaps[np.newaxis,:,:19]-fij_with_gaps[:,:,20,:19] # LxLx20 fj_i(a) = fj(b)-fij(-,b)
    fifj = fi_j[:,:,:,None]* fj_i[:,:,None,:]  # LxLx20x20 fifj(a,b) = (fi(a)-fij(a,-))*(fj(b)-fij(-,b))
    fifj_roll = np.transpose(fifj[:,:,:19,:19], axes=(0,2,1,3))
    fifj_reshape = fifj_roll.reshape((L*19,L*19))
    fij_roll = np.transpose(fij_with_gaps[:,:,:19,:19], axes=(0,2,1,3))
    fij_reshape = fij_roll.reshape((L*19,L*19))
    cov_matrix = gap_correction_matrix*(fij_reshape-gap_correction_matrix*fifj_reshape)
    return np.around(cov_matrix, np.finfo(np.double).precision)


def apply_norm_threshold(C, norm_threshold):
    L = C.shape[0]//19
    C_threshold = np.zeros_like(C)
    for i in range(L):
        for j in range(i,L):
            Cij = C[i*19:i*19+19,j*19:j*19+19]
            if np.linalg.norm(Cij)>=norm_threshold:
                C_threshold[i*19:i*19+19,j*19:j*19+19] = C[i*19:i*19+19,j*19:j*19+19]
                C_threshold[j*19:j*19+19,i*19:i*19+19] = np.transpose(C_threshold[i*19:i*19+19,j*19:j*19+19])
    return C_threshold


def get_norm_covariance_matrix_threshold(N):
    return 1/np.sqrt(0.85*N)


def covariance_matrix_for_uniform_distribution(L):
    U = np.zeros((L*19,L*19))
    for i in range(L):
        for a in range(19):
            for b in range(19):
                if (a==b):
                    U[i*19+a,i*19+b] = (1/20)*(1-1/20)
                else:
                    U[i*19+a,i*19+b] = -1/(20**2)
    return U

def shrink_covariance_matrix_towards_uniform(C, shrinkage_coeff):
    L = C.shape[0]//19
    U = covariance_matrix_for_uniform_distribution(L)
    return (1-shrinkage_coeff)*C + shrinkage_coeff*U

def prepare_covariance_matrix_for_zero_sum_gauge(C):
    L = C.shape[0]//19
    I = np.identity(C.shape[0])
    M = block_diag(*[np.ones((19,19))]*L)
    IM = I+M
    return np.matmul(IM, np.matmul(C, IM))
