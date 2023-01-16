import numpy as np
#import pkg_resources

def apply_uniform_pseudocounts_to_single_frequencies(fi, uniform_pc_rate):
    """ inputs single frequency in Lxq array, outputs single frequencies with uniform pseudocounts weighted by rate @uniform_pc_rate """
    q = fi.shape[1]
    return (1-uniform_pc_rate)*fi + uniform_pc_rate*np.full_like(fi, 1/q)


def apply_uniform_pseudocounts_to_covariance_matrix(C, fi_pc_without_gaps, tau_pc):
    """ inputs covariance matrix with shape (L*19,L*19) computed without gaps @C,
                single frequencies without gaps with pseudocounts already applied @fi_pc_without_gaps,
                desired pseudocount rate @tau_pc
        outputs covariance matrix computed with frequencies where uniform pseudocounts are applied """
    L = C.shape[0]//19
    C_pc = ((1-tau_pc)**2)*C
    for i in range(L):
        for a in range(19):
            for b in range(19):
                if (a!=b):
                    C_pc[i*19+a,i*19+b] = -fi_pc_without_gaps[i,a]*fi_pc_without_gaps[i,b]
                else:
                    C_pc[i*19+a,i*19+b] = fi_pc_without_gaps[i,a]*(1-fi_pc_without_gaps[i,a])
    return C_pc
