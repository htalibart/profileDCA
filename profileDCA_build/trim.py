import numpy as np

def get_positions_without_too_many_gaps(msa, max_gap_fraction_per_column, gap_symbol):
    """ returns a list of positions in @msa with fraction if @gap_symbol <= @max_gap_fraction_per_column """
    N, L = msa.shape
    nb_gaps_max = int(max_gap_fraction_per_column*N)
    keep_pos = []
    for i in range(L):
        if sum(msa[:,i]==gap_symbol)<=nb_gaps_max:
            keep_pos.append(i)
    return keep_pos

def trim_int_msa_array(msa, max_gap_fraction_per_column, gap_symbol):
    """
        @msa: int array
        @gap_symbol: int for gap symbol
        @max_gap_fraction_per_colum: max % gap per column in the MSA

        returns a smaller int array where all columns with > @max_gap_per_column are removed,
        and a list where l[index in smaller int array] = index in @msa
    """

    pos_to_keep = get_positions_without_too_many_gaps(msa, max_gap_fraction_per_column, gap_symbol)
    return msa[:,pos_to_keep], pos_to_keep


def insert_null_positions_to_complete_mrf_pos(mrf, sequence_length):
    new_mrf = {}
    new_mrf['mrf_pos_to_seq_pos'] = [pos for pos in range(sequence_length)]
    q = mrf['v'].shape[1]
    new_mrf['w'] = np.zeros((sequence_length,sequence_length,q,q))
    new_mrf['w'][np.ix_(list(mrf['mrf_pos_to_seq_pos']), list(mrf['mrf_pos_to_seq_pos']))] = mrf['w']
    new_mrf['v'] = np.zeros((sequence_length,q))
    new_mrf['v'][np.ix_(list(mrf['mrf_pos_to_seq_pos']))] = mrf['v']
    new_mrf['v_full'] = mrf['v']
    return new_mrf


