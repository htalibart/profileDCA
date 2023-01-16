import unittest
import pathlib, tempfile, shutil

from Bio import SeqIO

from profileDCA_build import __main__ as ppbuild_main
from profileDCA_utils import io_management as iom
from profileDCA_align.__main__ import *


def get_vi_conserved_letter(conserved_letter, via=1):
    via_conserved = via
    vi = np.ones(20)*(-via_conserved/(20-1))
    vi[conserved_letter] = via_conserved
    return vi

def get_fake_model(conserved_letters_list, ijabs=[], wijab=10, via=1):
    v = np.array([get_vi_conserved_letter(c, via=via) for c in conserved_letters_list])
    w = np.zeros((len(conserved_letters_list), len(conserved_letters_list), 20, 20))
    for ijab in ijabs:
        w[ijab] = wijab
        w[ijab[1],ijab[0],ijab[2],ijab[3]] = wijab
    return {'v':v, 'w':w}


class Test_PPalign(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_no_end_gap_free(self):
        potts_folders = [pathlib.Path(tempfile.mkdtemp()) for k in range(2)]
        with open(potts_folders[0]/"msa.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("AFC\n")
            of.write("> 0bis\n")
            of.write("ACC\n")
            for k in range(1,10):
                of.write('> '+str(k)+'\n')
                of.write("A--\n")
        with open(potts_folders[1]/"msa.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("AC\n")
            for k in range(1,10):
                of.write('> '+str(k)+'\n')
                of.write("A-\n")
        for pf in potts_folders:
            sequence_file = pf/"sequence.fasta"
            sequence = iom.get_first_sequence_in_fasta_file(pf/"msa.fasta")
            with open(sequence_file, 'w') as of:
                of.write("> 0\n")
                of.write(sequence)
        potts_objects = [ppbuild_main.build_mrf_to_folder(pf, pf/"sequence.fasta", pf/"msa.fasta", max_gap_v=1, max_gap_w=1, seq_id_threshold=1) for pf in potts_folders]
        output_folder = pathlib.Path(tempfile.mkdtemp())
        sequence_files = [pf/"sequence.fasta" for pf in potts_folders]
        align_objects_and_handle_files(potts_objects, output_folder, sequence_files, offset_v=1.25, no_free_end_gaps=True)
        output_fasta_file = output_folder/("aligned_sequences.fasta")
        aligned_records = list(SeqIO.parse(str(output_fasta_file), 'fasta'))
        aligned_sequences = [str(rec.seq) for rec in aligned_records]
        print(aligned_sequences)
        self.assertTrue(aligned_sequences==["AfC","A-C"])
        for pf in potts_folders:
            shutil.rmtree(pf)
        shutil.rmtree(output_folder)


#    def test_solver_convergence_problem(self): # DOES NOT CONVERGE!!
#        output_folder = pathlib.Path(tempfile.mkdtemp())
#        gap_open=1000
#        gap_extend=0
#        wijab=0.25
#        via=1
#        mrf1 = get_fake_model([0,1,2], ijabs=[(0,2,0,0)], wijab=wijab, via=via)
#        mrf2 = get_fake_model([0,2], ijabs=[(0,1,0,0)], wijab=wijab, via=via)
#        insert_opens = [np.array([0]*4), np.array([0,gap_open,0])]
#        insert_extends = [np.array([0]*4), np.array([0,gap_extend,0])]
#        mrfs = [mrf1, mrf2]
#        res = align_two_mrfs_with_solver(mrfs, output_folder, alpha_w=1, offset_v=0, gap_open=gap_open, gap_extend=gap_extend, insert_opens=insert_opens, insert_extends=insert_extends)


if __name__=='__main__':
    unittest.main()
