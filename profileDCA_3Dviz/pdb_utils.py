import numpy as np
from urllib.request import urlopen
from Bio import PDB

from profileDCA_utils import ppfunctions

d_3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def get_3to1(three):
    """ 3-letter amino acid code to 1-letter """
    if three in d_3to1:
        return d_3to1[three]
    else:
        raise Exception('Not an amino acid')

def fetch_pdb_file(pdb_id, outputfname):
    """ downloads .pdb file from PDB database """
    try:
        url = "https://files.rcsb.org/download/"+pdb_id+".pdb"
        pdbfile = urlopen(url)
        with open(str(outputfname)+".pdb",'wb') as output:
            output.write(pdbfile.read())
        return str(outputfname)+".pdb"
    except Exception as e:
        url = "https://files.rcsb.org/download/"+pdb_id+".cif"
        ciffile = urlopen(url)
        with open(str(outputfname)+".cif", 'wb') as output:
            output.write(ciffile.read())
        return str(outputfname)+".cif"

def get_pdb_chain(pdb_file, pdb_id, chain_id):
    """ opens PDB file @pdb_file with PDB id @pdb_id and outputs the BioPython chain object with chain id @chain_id """
    pdbfile = str(pdb_file)
    if pdbfile.endswith(".pdb"):
        structure = PDB.PDBParser().get_structure(pdb_id, pdbfile)
    elif pdbfile.endswith(".cif"):
        structure = Bio.PDB.MMCIFParser().get_structure(pdb_id, pdbfile)
    else:
        raise Exception("Unknown PDB file format")
    model = structure[0]
    pdb_chain = model[chain_id]
    return pdb_chain

def is_acceptable_residue(residue):
    """ inputs a BioPython residue object @residue, returns True if it is actually a residue """
    return residue.get_full_id()[3][0]==' '

def get_res_id(residue):
    """ inputs a BioPython residue object @residue, returns its id """
    return residue.get_full_id()[3][1]

def get_res_letter(residue):
    """ inputs a BioPython residue object @residue, returns its amino acid 1-letter code """
    return get_3to1(residue.get_resname())


def get_seq_pos_to_pdb_chain_pos(sequence, pdb_chain):
    """ inputs a protein sequence @sequence and a BioPython PDB chain object @pdb_chain
        returns mappings from sequence to PDB chain: @seq_to_pdb_chain where seq_to_pdb_chain[k] is the position in the PDB chain corresponding to position k in sequence """
    pdb_sequence_dict = {}
    for residue in pdb_chain:
        if is_acceptable_residue(residue):
            pdb_sequence_dict[get_res_id(residue)] = get_res_letter(residue)
    pdb_sequence = "".join(pdb_sequence_dict.values())
    seq_to_pdb_seq = ppfunctions.get_pos_first_seq_to_second_seq(sequence, pdb_sequence)
    seq_to_pdb_chain = []
    for seq_index in range(len(sequence)):
        pdb_seq_index = seq_to_pdb_seq[seq_index]
        if pdb_seq_index is None:
            pdb_chain_index = None
        else:
            pdb_chain_index = list(pdb_sequence_dict.keys())[pdb_seq_index]
            assert(get_res_letter(pdb_chain[pdb_chain_index])==sequence[seq_index])
        seq_to_pdb_chain.append(pdb_chain_index)
    return seq_to_pdb_chain


def get_mrf_pos_to_pdb_chain_pos(mrf_pos_to_seq_pos, sequence, pdb_chain):
    """ inputs a protein sequence @sequence, its mappings to a model @mrf_pos_to_seq_pos where mrf_pos_to_seq_pos[k] is the position in the model corresponding to position k in sequence, a BioPython PDB chain object @pdb_chain
        returns mappings from Potts model to PDB chain: @mrf_pos_to_pdb_chain_pos where mrf_pos_to_pdb_chain_pos[k] is the position in the PDB chain corresponding to position k in model """
    mrf_pos_to_pdb_chain_pos = []
    seq_pos_to_pdb_chain_pos = get_seq_pos_to_pdb_chain_pos(sequence, pdb_chain)
    for mrf_pos in range(len(mrf_pos_to_seq_pos)):
        seq_pos = mrf_pos_to_seq_pos[mrf_pos]
        if seq_pos is None:
            pdb_pos = None
        else:
            pdb_pos = seq_pos_to_pdb_chain_pos[seq_pos]
        mrf_pos_to_pdb_chain_pos.append(pdb_pos)
    return mrf_pos_to_pdb_chain_pos


def aa_distance(pos1, pos2, pdb_chain):
    """ returns euclidean distance between residues at positions @pos1 and @pos2 in structure represented by BioPython PDB chain object @pdb_chain """
    if (pos1 not in pdb_chain) or (pos2 not in pdb_chain):
        raise Exception("Position not in PDB chain")
    r1 = pdb_chain[pos1]
    r2 = pdb_chain[pos2]
    diff_vector = r1['CA'].coord - r2['CA'].coord
    return np.sqrt(np.sum(diff_vector*diff_vector))


def is_pdb_pair_contact(i_pdb, j_pdb, pdb_chain, contact_threshold=8):
    """ returns True if distance between residues at positions @i_pdb and @j_pdb in BioPython object @pdb_chain is lower than the threshold @contact_threshold in Angstrom """
    if (i_pdb is None) or (j_pdb is None):
        raise Exception("Position not in PDB chain")
    else:
        return aa_distance(i_pdb, j_pdb, pdb_chain)<=contact_threshold


def is_coupling_contact(i, j, pdb_chain, mrf_pos_to_pdb_chain_pos, contact_threshold=8):
    """ returns True if coupling between positions @i and @j in model is a True contact
        requires BioPython object @pdb_chain, its mappings to the model @mrf_pos_to_pdb_chain_pos and a threshold @contact_threshold to define a contact in Angstrom """
    i_pdb = mrf_pos_to_pdb_chain_pos[i]
    j_pdb = mrf_pos_to_pdb_chain_pos[j]
    return is_pdb_pair_contact(i_pdb, j_pdb, pdb_chain, contact_threshold=contact_threshold)
