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

    def test_align_without_missing_positions(self):
        potts_folders = [pathlib.Path(tempfile.mkdtemp()) for k in range(2)]
        with open(potts_folders[0]/"msa.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("ACDE\n")
            for k in range(2):
                of.write('> '+str(k)+'\n')
                of.write("ACDE\n")
            of.write("> other\n")
            of.write("FCDG")
        with open(potts_folders[1]/"msa.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("ADE\n")
            for k in range(2):
                of.write('> '+str(k)+'\n')
                of.write("ADE\n")
            of.write("> other\n")
            of.write("FDG")
        for pf in potts_folders:
            sequence_file = pf/"sequence.fasta"
            sequence = iom.get_first_sequence_in_fasta_file(pf/"msa.fasta")
            with open(sequence_file, 'w') as of:
                of.write("> 0\n")
                of.write(sequence)
        potts_objects = [ppbuild_main.build_mrf_to_folder(pf, pf/"sequence.fasta", pf/"msa.fasta", max_gap_v=0.1, max_gap_w=1, seq_id_threshold=1) for pf in potts_folders]
        output_folder = pathlib.Path(tempfile.mkdtemp())
        sequence_files = [pf/"sequence.fasta" for pf in potts_folders]
        align_objects_and_handle_files(potts_objects, output_folder, sequence_files)
        output_fasta_file = output_folder/("aligned_sequences.fasta")
        aligned_records = list(SeqIO.parse(str(output_fasta_file), 'fasta'))
        aligned_sequences = [str(rec.seq) for rec in aligned_records]
        print(aligned_sequences)
        assert(aligned_sequences==["AcDE","A-DE"])
        for pf in potts_folders:
            shutil.rmtree(pf)
        shutil.rmtree(output_folder)

    def test_realign_missing_positions(self):
        potts_folders = [pathlib.Path(tempfile.mkdtemp()) for k in range(2)]
        with open(potts_folders[0]/"msa.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("ACDE\n")
            for k in range(2):
                of.write('> '+str(k)+'\n')
                of.write("A--E\n")
            of.write("> other\n")
            of.write("F--G")
        with open(potts_folders[1]/"msa.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("ADE\n")
            for k in range(2):
                of.write('> '+str(k)+'\n')
                of.write("A-E\n")
            of.write("> other\n")
            of.write("F-G")
        for pf in potts_folders:
            sequence_file = pf/"sequence.fasta"
            sequence = iom.get_first_sequence_in_fasta_file(pf/"msa.fasta")
            with open(sequence_file, 'w') as of:
                of.write("> 0\n")
                of.write(sequence)
        potts_objects = [ppbuild_main.build_mrf_to_folder(pf, pf/"sequence.fasta", pf/"msa.fasta", max_gap_v=0.1, max_gap_w=1, seq_id_threshold=1) for pf in potts_folders]
        output_folder = pathlib.Path(tempfile.mkdtemp())
        sequence_files = [pf/"sequence.fasta" for pf in potts_folders]
        align_objects_and_handle_files(potts_objects, output_folder, sequence_files)
        output_fasta_file = output_folder/("aligned_sequences.fasta")
        aligned_records = list(SeqIO.parse(str(output_fasta_file), 'fasta'))
        aligned_sequences = [str(rec.seq) for rec in aligned_records]
        print(aligned_sequences)
        assert(aligned_sequences==["AcdE","A-dE"])
        for pf in potts_folders:
            shutil.rmtree(pf)
        shutil.rmtree(output_folder)


    def test_realign_missing_positions_2(self):
        potts_folders = [pathlib.Path(tempfile.mkdtemp()) for k in range(2)]
        with open(potts_folders[0]/"msa.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("ACDE\n")
            for k in range(2):
                of.write('> '+str(k)+'\n')
                of.write("A---\n")
            of.write("> other\n")
            of.write("F---")
        with open(potts_folders[1]/"msa.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("ADE\n")
            for k in range(2):
                of.write('> '+str(k)+'\n')
                of.write("A-E\n")
            of.write("> other\n")
            of.write("F-G")
        for pf in potts_folders:
            sequence_file = pf/"sequence.fasta"
            sequence = iom.get_first_sequence_in_fasta_file(pf/"msa.fasta")
            with open(sequence_file, 'w') as of:
                of.write("> 0\n")
                of.write(sequence)
        potts_objects = [ppbuild_main.build_mrf_to_folder(pf, pf/"sequence.fasta", pf/"msa.fasta", max_gap_v=0.1, max_gap_w=1, seq_id_threshold=1) for pf in potts_folders]
        output_folder = pathlib.Path(tempfile.mkdtemp())
        sequence_files = [pf/"sequence.fasta" for pf in potts_folders]
        align_objects_and_handle_files(potts_objects, output_folder, [pf/"sequence.fasta" for pf in potts_folders])
        output_fasta_file = output_folder/("aligned_sequences.fasta")
        aligned_records = list(SeqIO.parse(str(output_fasta_file), 'fasta'))
        aligned_sequences = [str(rec.seq) for rec in aligned_records]
        print(aligned_sequences)
        assert(aligned_sequences==["Acde","A-de"])
        for pf in potts_folders:
            shutil.rmtree(pf)
        shutil.rmtree(output_folder)


    def test_realign_missing_positions_3(self):
        potts_folders = [pathlib.Path(tempfile.mkdtemp()) for k in range(2)]
        with open(potts_folders[0]/"msa.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("ACD-\n")
            for k in range(2):
                of.write('> '+str(k)+'\n')
                of.write("A--E\n")
            of.write("> other\n")
            of.write("F--G")
        with open(potts_folders[0]/"sequence.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("ACD")
        with open(potts_folders[1]/"msa.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("ADW\n")
            for k in range(2):
                of.write('> '+str(k)+'\n')
                of.write("A-W\n")
            of.write("> other\n")
            of.write("F-G")
        with open(potts_folders[1]/"sequence.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("ADW")
        potts_objects = [ppbuild_main.build_mrf_to_folder(pf, pf/"sequence.fasta", pf/"msa.fasta", max_gap_v=0.1, max_gap_w=1, seq_id_threshold=1) for pf in potts_folders]
        output_folder = pathlib.Path(tempfile.mkdtemp())
        sequence_files = [pf/"sequence.fasta" for pf in potts_folders]
        align_objects_and_handle_files(potts_objects, output_folder, [pf/"sequence.fasta" for pf in potts_folders])
        output_fasta_file = output_folder/("aligned_sequences.fasta")
        aligned_records = list(SeqIO.parse(str(output_fasta_file), 'fasta'))
        aligned_sequences = [str(rec.seq) for rec in aligned_records]
        print(aligned_sequences)
        self.assertEqual(aligned_sequences,["Acd-","A-dw"])
        for pf in potts_folders:
            shutil.rmtree(pf)
        shutil.rmtree(output_folder)


    def test_zero_gap_open_at_trimmed(self):
        potts_folders = [pathlib.Path(tempfile.mkdtemp()) for k in range(2)]
        with open(potts_folders[0]/"msa.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("AWD\n")
            for k in range(4):
                of.write('> '+str(k)+'\n')
                of.write("A-D\n")
        with open(potts_folders[0]/"sequence.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("AWD")
        with open(potts_folders[1]/"msa.fasta", 'w') as of:
            for k in range(5):
                of.write('> '+str(k)+'\n')
                of.write("APWGHPD\n")
        with open(potts_folders[1]/"sequence.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("APWGHPD")
        potts_objects = [ppbuild_main.build_mrf_to_folder(pf, pf/"sequence.fasta", pf/"msa.fasta", max_gap_v=0.1, max_gap_w=1, seq_id_threshold=1) for pf in potts_folders]
        output_folder = pathlib.Path(tempfile.mkdtemp())
        sequence_files = [pf/"sequence.fasta" for pf in potts_folders]
        align_objects_and_handle_files(potts_objects, output_folder, [pf/"sequence.fasta" for pf in potts_folders])
        output_fasta_file = output_folder/("aligned_sequences.fasta")
        aligned_records = list(SeqIO.parse(str(output_fasta_file), 'fasta'))
        aligned_sequences = [str(rec.seq) for rec in aligned_records]
        print(aligned_sequences)
        self.assertEqual(aligned_sequences,["A-w---D","ApwghpD"])
        for pf in potts_folders:
            shutil.rmtree(pf)
        shutil.rmtree(output_folder)


    def test_zero_gap_open_at_trimmed_nothing_trimmed(self):
        potts_folders = [pathlib.Path(tempfile.mkdtemp()) for k in range(2)]
        with open(potts_folders[0]/"msa.fasta", 'w') as of:
            for k in range(5):
                of.write('> '+str(k)+'\n')
                of.write("ACD\n")
        with open(potts_folders[0]/"sequence.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("ACD")
        with open(potts_folders[1]/"msa.fasta", 'w') as of:
            for k in range(5):
                of.write('> '+str(k)+'\n')
                of.write("AECGHID\n")
        with open(potts_folders[1]/"sequence.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("AECGHID")
        potts_objects = [ppbuild_main.build_mrf_to_folder(pf, pf/"sequence.fasta", pf/"msa.fasta", max_gap_v=0.1, max_gap_w=1, seq_id_threshold=1) for pf in potts_folders]
        output_folder = pathlib.Path(tempfile.mkdtemp())
        gap_open=1000
        sequence_files = [pf/"sequence.fasta" for pf in potts_folders]
        align_objects_and_handle_files(potts_objects, output_folder, sequence_files, gap_open=gap_open)
        output_fasta_file = output_folder/("aligned_sequences.fasta")
        aligned_records = list(SeqIO.parse(str(output_fasta_file), 'fasta'))
        aligned_sequences = [str(rec.seq).upper() for rec in aligned_records]
        print(aligned_sequences)
        assert(aligned_sequences==["ACD----","AECGHID"] or aligned_sequences==["-ACD---","AECGHID"] or aligned_sequences==["----ACD","AECGHID"])
        for pf in potts_folders:
            shutil.rmtree(pf)
        shutil.rmtree(output_folder)

    def test_align_pos_not_same_trim(self):
        potts_folders = [pathlib.Path(tempfile.mkdtemp()) for k in range(2)]
        with open(potts_folders[0]/"msa.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("ACDE\n")
            for k in range(4):
                of.write('> '+str(k)+'\n')
                of.write("AC-E\n")
        with open(potts_folders[1]/"msa.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("ACDE\n")
            for k in range(2):
                of.write('> '+str(k)+'\n')
                of.write("A-DE\n")
        for pf in potts_folders:
            sequence_file = pf/"sequence.fasta"
            sequence = iom.get_first_sequence_in_fasta_file(pf/"msa.fasta")
            with open(sequence_file, 'w') as of:
                of.write("> 0\n")
                of.write(sequence)
        potts_objects = [ppbuild_main.build_mrf_to_folder(pf, pf/"sequence.fasta", pf/"msa.fasta", max_gap_v=0.1, max_gap_w=1, seq_id_threshold=1) for pf in potts_folders]
        output_folder = pathlib.Path(tempfile.mkdtemp())
        sequence_files = [pf/"sequence.fasta" for pf in potts_folders]
        align_objects_and_handle_files(potts_objects, output_folder, [pf/"sequence.fasta" for pf in potts_folders])
        output_fasta_file = output_folder/("aligned_sequences.fasta")
        aligned_records = list(SeqIO.parse(str(output_fasta_file), 'fasta'))
        aligned_sequences = [str(rec.seq) for rec in aligned_records]
        self.assertEqual(aligned_sequences,["AcdE","AcdE"])
        for pf in potts_folders:
            shutil.rmtree(pf)
        shutil.rmtree(output_folder)

    def test_align_pos_not_same_trim_2(self):
        potts_folders = [pathlib.Path(tempfile.mkdtemp()) for k in range(2)]
        with open(potts_folders[0]/"msa.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("ACDE\n")
            for k in range(4):
                of.write('> '+str(k)+'\n')
                of.write("AC-E\n")
        with open(potts_folders[1]/"msa.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("ACDFE\n")
            for k in range(2):
                of.write('> '+str(k)+'\n')
                of.write("A-DFE\n")
        for pf in potts_folders:
            sequence_file = pf/"sequence.fasta"
            sequence = iom.get_first_sequence_in_fasta_file(pf/"msa.fasta")
            with open(sequence_file, 'w') as of:
                of.write("> 0\n")
                of.write(sequence)
        potts_objects = [ppbuild_main.build_mrf_to_folder(pf, pf/"sequence.fasta", pf/"msa.fasta", max_gap_v=0.1, max_gap_w=1, seq_id_threshold=1) for pf in potts_folders]
        output_folder = pathlib.Path(tempfile.mkdtemp())
        sequence_files = [pf/"sequence.fasta" for pf in potts_folders]
        align_objects_and_handle_files(potts_objects, output_folder, [pf/"sequence.fasta" for pf in potts_folders])
        output_fasta_file = output_folder/("aligned_sequences.fasta")
        aligned_records = list(SeqIO.parse(str(output_fasta_file), 'fasta'))
        aligned_sequences = [str(rec.seq) for rec in aligned_records]
        self.assertEqual(aligned_sequences,["Acd-E","AcdfE"])
        for pf in potts_folders:
            shutil.rmtree(pf)
        shutil.rmtree(output_folder)

    def test_align_conserved_insertion_and_small_gap_open(self):
        potts_folders = [pathlib.Path(tempfile.mkdtemp()) for k in range(2)]
        with open(potts_folders[0]/"msa.fasta", 'w') as of:
            of.write('> 0_\n')
            of.write("AC\n")
            for k in range(4):
                of.write('> '+str(k)+'\n')
                of.write("AC\n")
        with open(potts_folders[1]/"msa.fasta", 'w') as of:
            of.write('> 0_\n')
            of.write("ADC\n")
            for k in range(2):
                of.write('> '+str(k)+'\n')
                of.write("AD-\n")
        for pf in potts_folders:
            sequence_file = pf/"sequence.fasta"
            sequence = iom.get_first_sequence_in_fasta_file(pf/"msa.fasta")
            with open(sequence_file, 'w') as of:
                of.write("> 0\n")
                of.write(sequence)
        potts_objects = [ppbuild_main.build_mrf_to_folder(pf, pf/"sequence.fasta", pf/"msa.fasta", max_gap_v=0.1, max_gap_w=1, seq_id_threshold=1) for pf in potts_folders]
        output_folder = pathlib.Path(tempfile.mkdtemp())
        sequence_files = [pf/"sequence.fasta" for pf in potts_folders]
        align_objects_and_handle_files(potts_objects, output_folder, sequence_files, offset_v=1.25, gap_open=1)
        output_fasta_file = output_folder/("aligned_sequences.fasta")
        aligned_records = list(SeqIO.parse(str(output_fasta_file), 'fasta'))
        aligned_sequences = [str(rec.seq) for rec in aligned_records]
        print(aligned_sequences)
        self.assertEqual(aligned_sequences,["A-c","Adc"])
        for pf in potts_folders:
            shutil.rmtree(pf)
        shutil.rmtree(output_folder)


    def test_begin_gap_free_for_different_letters(self):
        potts_folders = [pathlib.Path(tempfile.mkdtemp()) for k in range(2)]
        with open(potts_folders[0]/"msa.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("AFCE\n")
            for k in range(1,10):
                of.write('> '+str(k)+'\n')
                of.write("---E\n")
        with open(potts_folders[1]/"msa.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("DE\n")
            for k in range(1,10):
                of.write('> '+str(k)+'\n')
                of.write("-E\n")
        for pf in potts_folders:
            sequence_file = pf/"sequence.fasta"
            sequence = iom.get_first_sequence_in_fasta_file(pf/"msa.fasta")
            with open(sequence_file, 'w') as of:
                of.write("> 0\n")
                of.write(sequence)
        potts_objects = [ppbuild_main.build_mrf_to_folder(pf, pf/"sequence.fasta", pf/"msa.fasta", max_gap_v=0.1, max_gap_w=1, seq_id_threshold=1) for pf in potts_folders]
        output_folder = pathlib.Path(tempfile.mkdtemp())
        sequence_files = [pf/"sequence.fasta" for pf in potts_folders]
        align_objects_and_handle_files(potts_objects, output_folder, sequence_files, offset_v=1.25, gap_open=1)
        output_fasta_file = output_folder/("aligned_sequences.fasta")
        aligned_records = list(SeqIO.parse(str(output_fasta_file), 'fasta'))
        aligned_sequences = [str(rec.seq) for rec in aligned_records]
        print(aligned_sequences)
        self.assertEqual(aligned_sequences,["afc-E","---dE"])
        for pf in potts_folders:
            shutil.rmtree(pf)
        shutil.rmtree(output_folder)


    def test_begin_gap_free_for_similar_letters(self):
        potts_folders = [pathlib.Path(tempfile.mkdtemp()) for k in range(2)]
        with open(potts_folders[0]/"msa.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("AFLE\n")
            for k in range(1,10):
                of.write('> '+str(k)+'\n')
                of.write("---E\n")
        with open(potts_folders[1]/"msa.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("VE\n")
            for k in range(1,10):
                of.write('> '+str(k)+'\n')
                of.write("-E\n")
        for pf in potts_folders:
            sequence_file = pf/"sequence.fasta"
            sequence = iom.get_first_sequence_in_fasta_file(pf/"msa.fasta")
            with open(sequence_file, 'w') as of:
                of.write("> 0\n")
                of.write(sequence)
        potts_objects = [ppbuild_main.build_mrf_to_folder(pf, pf/"sequence.fasta", pf/"msa.fasta", max_gap_v=0.1, max_gap_w=1, seq_id_threshold=1) for pf in potts_folders]
        output_folder = pathlib.Path(tempfile.mkdtemp())
        sequence_files = [pf/"sequence.fasta" for pf in potts_folders]
        align_objects_and_handle_files(potts_objects, output_folder, sequence_files, offset_v=1.25)
        output_fasta_file = output_folder/("aligned_sequences.fasta")
        aligned_records = list(SeqIO.parse(str(output_fasta_file), 'fasta'))
        aligned_sequences = [str(rec.seq) for rec in aligned_records]
        print(aligned_sequences)
        self.assertEqual(aligned_sequences,["aflE","--vE"])
        for pf in potts_folders:
            shutil.rmtree(pf)
        shutil.rmtree(output_folder)


    def test_end_gap_free(self):
        potts_folders = [pathlib.Path(tempfile.mkdtemp()) for k in range(2)]
        with open(potts_folders[0]/"msa.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("AFCE\n")
            for k in range(1,10):
                of.write('> '+str(k)+'\n')
                of.write("A---\n")
        with open(potts_folders[1]/"msa.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("AD\n")
            for k in range(1,10):
                of.write('> '+str(k)+'\n')
                of.write("A-\n")
        for pf in potts_folders:
            sequence_file = pf/"sequence.fasta"
            sequence = iom.get_first_sequence_in_fasta_file(pf/"msa.fasta")
            with open(sequence_file, 'w') as of:
                of.write("> 0\n")
                of.write(sequence)
        potts_objects = [ppbuild_main.build_mrf_to_folder(pf, pf/"sequence.fasta", pf/"msa.fasta", max_gap_v=0.1, max_gap_w=1, seq_id_threshold=1) for pf in potts_folders]
        output_folder = pathlib.Path(tempfile.mkdtemp())
        sequence_files = [pf/"sequence.fasta" for pf in potts_folders]
        align_objects_and_handle_files(potts_objects, output_folder, sequence_files, offset_v=1.25)
        output_fasta_file = output_folder/("aligned_sequences.fasta")
        aligned_records = list(SeqIO.parse(str(output_fasta_file), 'fasta'))
        aligned_sequences = [str(rec.seq) for rec in aligned_records]
        print(aligned_sequences)
        self.assertEqual(aligned_sequences,["Afce","Ad--"])
        for pf in potts_folders:
            shutil.rmtree(pf)
        shutil.rmtree(output_folder)


    def test_begin_gap_free_with_unaligned_column(self):
        potts_folders = [pathlib.Path(tempfile.mkdtemp()) for k in range(2)]
        with open(potts_folders[0]/"msa.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("AFCE\n")
            for k in range(1,10):
                of.write('> '+str(k)+'\n')
                of.write("---E\n")
        with open(potts_folders[1]/"msa.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("GDE\n")
            for k in range(1,10):
                of.write('> '+str(k)+'\n')
                of.write("G-E\n")
        for pf in potts_folders:
            sequence_file = pf/"sequence.fasta"
            sequence = iom.get_first_sequence_in_fasta_file(pf/"msa.fasta")
            with open(sequence_file, 'w') as of:
                of.write("> 0\n")
                of.write(sequence)
        potts_objects = [ppbuild_main.build_mrf_to_folder(pf, pf/"sequence.fasta", pf/"msa.fasta", max_gap_v=0.1, max_gap_w=1, seq_id_threshold=1) for pf in potts_folders]
        output_folder = pathlib.Path(tempfile.mkdtemp())
        sequence_files = [pf/"sequence.fasta" for pf in potts_folders]
        align_objects_and_handle_files(potts_objects, output_folder, sequence_files, offset_v=1.25)
        output_fasta_file = output_folder/("aligned_sequences.fasta")
        aligned_records = list(SeqIO.parse(str(output_fasta_file), 'fasta'))
        aligned_sequences = [str(rec.seq) for rec in aligned_records]
        print(aligned_sequences)
        self.assertEqual(aligned_sequences,["afcE","-gdE"])
        for pf in potts_folders:
            shutil.rmtree(pf)
        shutil.rmtree(output_folder)


    def test_end_gap_free_with_unaligned_column(self):
        potts_folders = [pathlib.Path(tempfile.mkdtemp()) for k in range(2)]
        with open(potts_folders[0]/"msa.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("ADCE\n")
            for k in range(1,10):
                of.write('> '+str(k)+'\n')
                of.write("-DC-\n")
        with open(potts_folders[1]/"msa.fasta", 'w') as of:
            of.write('> 0\n')
            of.write("DE\n")
            for k in range(1,10):
                of.write('> '+str(k)+'\n')
                of.write("DE\n")
        for pf in potts_folders:
            sequence_file = pf/"sequence.fasta"
            sequence = iom.get_first_sequence_in_fasta_file(pf/"msa.fasta")
            with open(sequence_file, 'w') as of:
                of.write("> 0\n")
                of.write(sequence)
        potts_objects = [ppbuild_main.build_mrf_to_folder(pf, pf/"sequence.fasta", pf/"msa.fasta", max_gap_v=0.1, max_gap_w=1, seq_id_threshold=1) for pf in potts_folders]
        output_folder = pathlib.Path(tempfile.mkdtemp())
        sequence_files = [pf/"sequence.fasta" for pf in potts_folders]
        align_objects_and_handle_files(potts_objects, output_folder, sequence_files, offset_v=1.25)
        output_fasta_file = output_folder/("aligned_sequences.fasta")
        aligned_records = list(SeqIO.parse(str(output_fasta_file), 'fasta'))
        aligned_sequences = [str(rec.seq) for rec in aligned_records]
        print(aligned_sequences)
        self.assertEqual(aligned_sequences,["aDCe","-DE-"])
        for pf in potts_folders:
            shutil.rmtree(pf)
        shutil.rmtree(output_folder)

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
        potts_objects = [ppbuild_main.build_mrf_to_folder(pf, pf/"sequence.fasta", pf/"msa.fasta", max_gap_v=1, max_gap_w=0, seq_id_threshold=1) for pf in potts_folders]
        output_folder = pathlib.Path(tempfile.mkdtemp())
        sequence_files = [pf/"sequence.fasta" for pf in potts_folders]
        align_objects_and_handle_files(potts_objects, output_folder, sequence_files, offset_v=1.25, no_free_end_gaps=True)
        output_fasta_file = output_folder/("aligned_sequences.fasta")
        aligned_records = list(SeqIO.parse(str(output_fasta_file), 'fasta'))
        aligned_sequences = [str(rec.seq) for rec in aligned_records]
        print(aligned_sequences)
        self.assertTrue(aligned_sequences==["AfC","A-C"])
        align_objects_and_handle_files(potts_objects, output_folder, sequence_files, offset_v=1.25, no_free_end_gaps=False)
        output_fasta_file = output_folder/("aligned_sequences.fasta")
        aligned_records = list(SeqIO.parse(str(output_fasta_file), 'fasta'))
        aligned_sequences = [str(rec.seq) for rec in aligned_records]
        print(aligned_sequences)
        self.assertEqual(aligned_sequences,["AFc","AC-"])

        for pf in potts_folders:
            shutil.rmtree(pf)
        shutil.rmtree(output_folder)

if __name__=='__main__':
    unittest.main()
