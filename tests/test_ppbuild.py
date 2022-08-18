import unittest
import pathlib, tempfile, shutil
import numpy as np

from profileDCA_build import msa_statistics, mrf_inference, covariance_processing, pseudocounts
from profileDCA_build import __main__ as ppbuild_main
from profileDCA_utils import io_management as iom

from profileDCA_viz import ppviz

HHBLITS_DATABASE = "/home/htalibart/databases/UniRef30_2020_06_hhsuite/UniRef30_2020_06"

class Test_PPbuild(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_covariance_matrix_gap_correction_conserved(self):
        msa_file_with_gaps = pathlib.Path("/tmp/")/(next(tempfile._get_candidate_names())+'.fasta')
        msa_file_without_gaps = pathlib.Path("/tmp/")/(next(tempfile._get_candidate_names())+'.fasta')
        sequences_without_gaps = ['AA']*6
        sequences_with_gaps = sequences_without_gaps+['A-']*2+['--']*2
        with open(msa_file_with_gaps, 'w') as of:
            for index, seq in enumerate(sequences_with_gaps):
                of.write('>'+str(index)+'\n')
                of.write(seq+'\n')
        with open(msa_file_without_gaps, 'w') as of:
            for index, seq in enumerate(sequences_without_gaps):
                of.write('>'+str(index)+'\n')
                of.write(seq+'\n')
        msa = iom.get_int_msa_array(msa_file_without_gaps)
        fi, fij = msa_statistics.compute_frequencies(msa)
        C = covariance_processing.compute_covariance_matrix_without_gaps(fi, fij)
        self.assertTrue(not np.any(C))
        msa = iom.get_int_msa_array(msa_file_with_gaps)
        fi, fij = msa_statistics.compute_frequencies(msa)
        C = covariance_processing.compute_covariance_matrix_without_gaps(fi, fij)
        self.assertTrue(np.allclose(C, np.zeros_like(C)))
        msa_file_with_gaps.unlink()
        msa_file_without_gaps.unlink()

    def test_covariance_matrix_gap_correction_coupling(self):
        msa_file_with_gaps = pathlib.Path("/tmp/")/(next(tempfile._get_candidate_names())+'.fasta')
        msa_file_without_gaps = pathlib.Path("/tmp/")/(next(tempfile._get_candidate_names())+'.fasta')
        sequences_without_gaps = ['AR']*3+['ND']*3
        sequences_with_gaps = sequences_without_gaps+['A-']*2+['--']*2
        with open(msa_file_with_gaps, 'w') as of:
            for index, seq in enumerate(sequences_with_gaps):
                of.write('>'+str(index)+'\n')
                of.write(seq+'\n')
        with open(msa_file_without_gaps, 'w') as of:
            for index, seq in enumerate(sequences_without_gaps):
                of.write('>'+str(index)+'\n')
                of.write(seq+'\n')
        msa = iom.get_int_msa_array(msa_file_without_gaps)
        fi, fij = msa_statistics.compute_frequencies(msa)
        C_nogap = covariance_processing.compute_covariance_matrix_without_gaps(fi, fij)
        self.assertTrue(np.any(C_nogap))
        msa = iom.get_int_msa_array(msa_file_with_gaps)
        fi, fij = msa_statistics.compute_frequencies(msa)
        C = covariance_processing.compute_covariance_matrix_without_gaps(fi, fij)
        fi_without_gaps = msa_statistics.fi_with_gaps_to_fi_without_gaps(fi)
        # TEST PSEUDO-COUNTS
        uniform_pc_rate=0
        fi_pc_without_gaps = pseudocounts.apply_uniform_pseudocounts_to_single_frequencies(fi_without_gaps, uniform_pc_rate)
        C_pc = pseudocounts.apply_uniform_pseudocounts_to_covariance_matrix(C, fi_pc_without_gaps, uniform_pc_rate)
        self.assertTrue(np.array_equal(C, C_pc))
        uniform_pc_rate=0.5
        fi_pc_without_gaps = pseudocounts.apply_uniform_pseudocounts_to_single_frequencies(fi_without_gaps, uniform_pc_rate)
        C_pc = pseudocounts.apply_uniform_pseudocounts_to_covariance_matrix(C, fi_pc_without_gaps, uniform_pc_rate)
        self.assertFalse(np.array_equal(C, C_pc))
        msa_file_with_gaps.unlink()
        msa_file_without_gaps.unlink()


    def test_covariance_matrix_threshold(self):
        L = 10
        C = np.ones((L*19,L*19))
        C[0*19+0,1*19+0]=0.5
        C[1*19+0,0*19+0]=0.5
        self.assertTrue(np.array_equal(C, C.T))
        self.assertTrue(np.array_equal(covariance_processing.apply_norm_threshold(C,0),C))
        self.assertTrue(np.array_equal(covariance_processing.apply_norm_threshold(C,1000),np.zeros_like(C)))

    def test_build_model_from_msa(self):
        files_dict = {}
        files_dict['sequence'] = pathlib.Path("/tmp/")/(next(tempfile._get_candidate_names())+'.fasta')
        sequence = "AEIKHYQFNVVMTCSGCSGAVNKVLTKLEPDVSKIDISLEKQLVDVYTTLPYDFILEKIKKTGKEVRSGKQL"
        with open(files_dict['sequence'], 'w') as sf:
            sf.write('> 1CC8\n')
            sf.write(sequence+"\n")
        files_dict['msa'] = pathlib.Path("/tmp/")/(next(tempfile._get_candidate_names())+'.fasta')
        with open(files_dict['msa'], 'w') as sf:
            sf.write('> 1CC8\n')
            sf.write('AEIKHYQFNVVMTCSGCSGAVNKVLTKLEPDVSKIDISLEKQLVDVYTTLPYDFILEKIKKTGKEVRSGKQL\n')
            sf.write('> 1CC8_2\n')
            sf.write('-EIKHYQFNVVMTCSGCSGAVNKVLTKLEPDVSKIDISLEKQLVDVYTTLPYDFILEKI--TGKEVRSGKQL\n')
        #output_file = pathlib.Path("/tmp/")/(next(tempfile._get_candidate_names())+'.fasta')
        #seq_id_threshold = 80
        potts_folder = pathlib.Path(tempfile.mkdtemp())
        mrf = ppbuild_main.get_mrf_from_processed_msa(files_dict['msa'], max_gap_v=0.4)
        self.assertEqual(mrf['v'].shape[1],mrf['v_full'].shape[1])
        self.assertEqual(mrf['v_full'].shape[0],len(sequence))
        self.assertEqual(mrf['v'].shape[0],len(sequence)-3)

        iom.mrf_to_files(mrf, potts_folder)
        mrf_from_folder = iom.mrf_from_folder(potts_folder)

        for key in files_dict:
            files_dict[key].unlink()

        shutil.rmtree(potts_folder)


    def test_build_model_from_sequence(self):
        files_dict = {}
        files_dict['sequence'] = pathlib.Path("/tmp/")/(next(tempfile._get_candidate_names())+'.fasta')
        sequence = "AEIKHYQFNVVMTCSGCSGAVNKVLTKLEPDVSKIDISLEKQLVDVYTTLPYDFILEKIKKTGKEVRSGKQL"
        with open(files_dict['sequence'], 'w') as sf:
            sf.write('> 1CC8\n')
            sf.write(sequence+"\n")
        potts_folder = pathlib.Path(tempfile.mkdtemp())
        mrf = ppbuild_main.build_mrf_to_folder(potts_folder, files_dict['sequence'], msa_file=None, database=HHBLITS_DATABASE)
        #mrf = iom.mrf_from_folder(potts_folder)
        ppviz.visualize_mrf(mrf)
        shutil.rmtree(potts_folder)



if __name__=='__main__':
    unittest.main()
