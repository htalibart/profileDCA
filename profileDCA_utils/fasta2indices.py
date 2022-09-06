#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import sys
import csv
import pathlib
import numpy as np
from Bio import AlignIO
ExtendedIUPACProtein='ACDEFGHIKLMNPQRSTVWYBXZJUO'

from profileDCA_utils import manage_positions
from profileDCA_utils import io_management as iom

def fasta2indices(alignment_handle, start_at_1=False):
    align = AlignIO.read(alignment_handle, 'fasta')
    tuple_list=[]
    assert(len(align) == 2)
    if start_at_1:
        pos0, pos1 = 0, 0
    else:
        pos0, pos1 = -1, -1
    for pair in zip(align[0].seq, align[1].seq):
        aligned_residues = True
        if pair[0].upper() in ExtendedIUPACProtein:
            pos0 += 1
        else:
            aligned_residues = False
        if pair[1].upper() in ExtendedIUPACProtein:
            pos1 += 1
        else:
            aligned_residues = False
        if aligned_residues:
            tuple_list.append((pos0, pos1))
    return tuple_list


def fasta2csv(alignment_handle, output_handle, seq_pos_to_mrf_pos=None, start_at_1=False):
    tuple_list = fasta2indices(alignment_handle, start_at_1=start_at_1)
    csv_writer = csv.writer(output_handle)
    csv_writer.writerow(['pos_ref','pos_2'])
    for pair in tuple_list:
        if seq_pos_to_mrf_pos is not None:
            pair = [seq_pos_to_mrf_pos[k][pair[k]] for k in range(2)]
        if (pair[0] is not None) and (pair[1] is not None):
            csv_writer.writerow(pair)
        else:
            print(pair, "was not aligned")

def fasta2indices_with_gaps(alignment_handle, start_at_1=False):
    align = AlignIO.read(alignment_handle, 'fasta')
    tuple_list=[]
    assert(len(align) == 2)
    if start_at_1:
        pos0, pos1 = 0, 0
    else:
        pos0, pos1 = -1, -1
    for pair in zip(align[0].seq, align[1].seq):
        if pair[0].upper() in ExtendedIUPACProtein:
            pos0+=1
            if pair[1].upper() in ExtendedIUPACProtein:
                pos1+=1
                tuple_list.append((pos0,pos1))
            else:
                tuple_list.append((pos0,'-'))
        else:
            if pair[1].upper() in ExtendedIUPACProtein:
                pos1+=1
                tuple_list.append(('-',pos1))
            else:
                tuple_list.append(('-','-'))
    return tuple_list


def fasta2csv_with_gaps(alignment_handle, output_handle, seq_pos_to_mrf_pos=None, start_at_1=False):
    tuple_list = fasta2indices_with_gaps(alignment_handle, start_at_1=start_at_1)
    csv_writer = csv.writer(output_handle)
    csv_writer.writerow(['pos_ref','pos_2'])
    for pair in tuple_list:
        if seq_pos_to_mrf_pos is not None:
            pair = []
            for k in range(2):
                pos_in_mrf = pair[k]
                if pos_in_mrf=='-':
                    pos_in_seq='-'
                else:
                    pos_in_seq = seq_pos_to_mrf_pos[k][pos_in_mrf]
                pair.append(pos_in_seq)
        if (pair[0] is not None) and (pair[1] is not None):
            csv_writer.writerow(pair)
        else:
            print(pair, "was not aligned")


def main(args=sys.argv[1:]):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'input_alignment', help='Pairwise alignment in fasta format',
        type=argparse.FileType('r')
    )
    parser.add_argument(
        'output_csv', help='Output file in (csv)',
        type=argparse.FileType('w')
    )
    parser.add_argument('-pf', '--potts_folders', help="Potts folders", type=pathlib.Path, nargs='+', default=[])
    args = vars(parser.parse_args(args))

    
    seq_pos_to_mrf_pos_list = [manage_positions.reverse_pos1_to_pos2(np.load(potts_folder/"mrf_pos_to_seq_pos.npy"), len(iom.get_first_sequence_in_fasta_file(potts_folder/"sequence.fasta"))) for potts_folder in args["potts_folders"]]
    fasta2csv(args['input_alignment'], args["output_csv"], seq_pos_to_mrf_pos_list)


if __name__ == '__main__':
    main()
