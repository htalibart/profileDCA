from math import *
import numpy as np


def simulate_uniform_pc_on_v(v, v_rescaling_tau=0.5, **kwargs):
    q = v.shape[1]
    resc_v = np.zeros_like(v)
    for i in range(len(v)):
        vi = v[i]
        resc_vi = np.zeros_like(vi)
        S = 0
        for b in range(q):
            S+=exp(vi[b])

        resc_tmp= np.zeros_like(vi)
        for a in range(q):
            resc_tmp[a] = log((1-v_rescaling_tau)*exp(vi[a])/S + v_rescaling_tau/q)

        resc_vi = np.zeros_like(vi)
        S_all = np.sum(resc_tmp)
        for a in range(q):
            resc_vi[a] = resc_tmp[a]-(1/q)*S_all
        resc_v[i] = resc_vi
    return resc_v



def simulate_uniform_pc_on_wij(w, rescaling_tau=0.5, beta=10, **kwargs):
    w_flat = w.flatten()
    S = sum([exp(beta*elt) for elt in w_flat])
    
    resc_tmp = np.zeros_like(w_flat)
    for elt in range(len(w_flat)):
        frac = exp(beta*w_flat[elt])/S
        resc_tmp[elt] = log((1-rescaling_tau)*frac+rescaling_tau/len(w_flat))
    
    resc_flat = np.zeros_like(w_flat)
    S_all = np.sum(resc_tmp)
    for elt in range(len(w_flat)):
        resc_flat[elt] = (1/beta)*(resc_tmp[elt]-(1/len(w_flat))*S_all)

    resc_unflat = resc_flat.reshape(w.shape)
    
    return resc_unflat




def simulate_uniform_pc_on_w(w, w_rescaling_tau=0.5, beta_softmax_w=10, **kwargs):
    resc_w = np.zeros_like(w)
    for i in range(len(resc_w)-1):
        for j in range(i+1,len(resc_w)):
            resc_w[i,j] = simulate_uniform_pc_on_wij(w[i][j], rescaling_tau=w_rescaling_tau, beta=beta_softmax_w)
            resc_w[j,i]=np.transpose(resc_w[i,j])
    return resc_w
