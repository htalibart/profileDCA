import unittest
import pathlib, tempfile, shutil

from Bio import SeqIO

from profileDCA_build import __main__ as ppbuild_main
from profileDCA_utils import io_management as iom
from profileDCA_align.__main__ import *

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

if __name__=='__main__':
    unittest.main()
