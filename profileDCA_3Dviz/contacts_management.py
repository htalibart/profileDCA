import tempfile
import pathlib
import numpy as np
from collections import OrderedDict
from itertools import islice
import matplotlib.pyplot as plt

from profileDCA_3Dviz import top_couplings, pdb_utils
from profileDCA_utils import ppfunctions

def get_contact_scores_for_w(w):
    """ returns an ordered dictionary of contact scores for positions in the train MSA """
    w_norms = ppfunctions.compute_w_norms(w)
    top_ind1, top_ind2 = top_couplings.get_top_pairs(w_norms, reverse_order=False)
    contact_scores = OrderedDict()
    for ind1, ind2 in zip(top_ind1, top_ind2):
        contact_scores[(ind1, ind2)] = w_norms[ind1, ind2]
    return contact_scores


def get_contact_scores_with_pdb_indexes(mrf, original_sequence, pdb_chain):
    """ returns an ordered dictionary of contact scores for positions in the PDB chain """
    mrf_pos_to_pdb_pos = pdb_utils.get_mrf_pos_to_pdb_chain_pos(mrf['mrf_pos_to_seq_pos'], original_sequence, pdb_chain)
    pm_scores = get_contact_scores_for_w(mrf['w'])
    pdb_contact_scores = OrderedDict()
    for c in pm_scores:
        c_pdb = (mrf_pos_to_pdb_pos[c[0]], mrf_pos_to_pdb_pos[c[1]])
        if (c_pdb[0] is not None) and (c_pdb[1] is not None):
            pdb_contact_scores[c_pdb] = pm_scores[c]
    return pdb_contact_scores


def get_colored_true_false_dicts(pdb_couplings_dict, pdb_chain, colors={True:'blue', False:'red'}, contact_threshold=8):
    """ tf_d['blue'][c] is the strength of coupling c predicted as true contact and tf_d['red'][c] is the strength of coupling c not predicted as contact """
    tf_d = {colors[val]:OrderedDict() for val in colors}
    for pdb_c in pdb_couplings_dict:
        if not None in pdb_c:
            tf_d[colors[is_true_contact(pdb_c, pdb_chain, contact_threshold=contact_threshold)]][pdb_c] = pdb_couplings_dict[pdb_c]
    return tf_d

def remove_couplings_too_close(couplings_dict, coupling_sep_min):
    """ returns dictionary where couplings separated by less than @coupling_sep_min are removed """
    ok_dict = OrderedDict()
    for c in couplings_dict:
        if None not in c:
            if abs(c[0]-c[1])>coupling_sep_min:
                ok_dict[c] = couplings_dict[c]
    return ok_dict

def get_smaller_dict(couplings_dict, nb_couplings):
    """ returns couplings dictionary with only the strongest @nb_couplings couplings """
    new_dict = OrderedDict()
    for c in couplings_dict:
        if len(new_dict)<nb_couplings:
            new_dict[c] = couplings_dict[c]
    return new_dict


def get_cutoff_smaller_than(couplings_dict, score_cutoff):
    """ returns index at which coupling strength starts to be lower than @score_cutoff in ordered @couplings_dict """
    y = list(couplings_dict.values())
    if y[0]<score_cutoff:
        return 0
    else:
        ind=0
        while (y[ind]>=score_cutoff):
            ind+=1
        return ind

def get_exclus_overlaps(couplings_dict_, tops):
    """ returns a list of 3 dicts : [couplings_dict[0]\couplings_dict[1], couplings_dict[0] inter couplings_dict[1], couplings_dict[1]\couplings_dict[0]. Overlap -> mean score """
    if len(couplings_dict_)==1:
        return [OrderedDict(islice(couplings_dict_[0].items(), 0, tops[0]))]
    else:
        couplings_dict = [OrderedDict(islice(couplings_dict_[k].items(), 0, int(tops[k]))) for k in range(2)]
        exclus = []
        overlaps = {}
        for k in range(2):
            exclus.append(OrderedDict())
            for c in couplings_dict[k]:
                if not ( (c in couplings_dict[(k+1)%2]) or ((c[1],c[0]) in couplings_dict[(k+1)%2]) ):
                    exclus[k][c] = couplings_dict[k][c]
                elif c in couplings_dict[(k+1)%2]:
                    overlaps[frozenset(c)] =  (couplings_dict[k][c]+couplings_dict[(k+1)%2][c])/2
                elif (c[1],c[0]) in couplings_dict[(k+1)%2]:
                    overlaps[frozenset(c)] = (couplings_dict[k][c]+couplings_dict[(k+1)%2][(c[1],c[0])])/2
        overlaps_ordered = OrderedDict({tuple(k): v for k, v in sorted(overlaps.items(), key=lambda item: item)}) # order by mean score
    return exclus+[overlaps_ordered]


def get_normalized_ordered_dict(od):
    """ normalizes dictionary values between 0 and 1 """
    fact=1/sum(od.values())
    new_od = OrderedDict()
    for key in od:
        new_od[key] = od[key]*fact
    return new_od
