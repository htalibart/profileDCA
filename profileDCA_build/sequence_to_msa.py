import pathlib, tempfile, os, subprocess
from profileDCA_build import msa_processing

def get_a3m_from_hhblits(input_file, output_file, database, maxfilt=100000, realign_max=100000, B=100000, Z=100000, n=3, e=0.001, retry_hhblits_with_memory_limit_if_fail=False, hhr_file=None):
    """ calls HH-blits with arguments recommended for CCMpred : https://github.com/soedinglab/CCMpred/wiki/FAQ """
    if hhr_file is None:
        hhr_file = pathlib.Path('.'.join(str(output_file).split('.')[:-1])+".hhr")
    hhblits_call = "hhblits -maxfilt "+str(maxfilt)+" -realign_max "+str(realign_max)+" -d "+str(database)+" -all -B "+str(B)+" -Z "+str(Z)+" -n "+str(n)+" -e "+str(e)+" -i "+str(input_file)+" -oa3m "+str(output_file)+" -o "+str(hhr_file)
    print(hhblits_call)
    subprocess.Popen(hhblits_call, shell=True).wait()
    if not output_file.exists() and retry_hhblits_with_memory_limit_if_fail:
        print("HHblits failed for some reason, trying again with a memory limit")
        memory_friendly_call = hhblits_call+" -cpu 1 -maxmem 1"
        subprocess.Popen(memory_friendly_call, shell=True).wait()
        if not output_file.exists():
            raise Exception("HHblits failed. Protein is probably too long ?")
    return hhr_file


def call_reformat(input_file, output_file):
    call = "reformat.pl a3m fas "+str(input_file)+" "+str(output_file)+" -r"
    print(call)
    subprocess.Popen(call, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).wait()
    if not output_file.is_file():
        raise Exception("Reformat failed")

def get_fasta_from_hhblits(input_file, output_file, database, maxfilt=100000, realign_max=100000, B=100000, Z=100000, n=3, e=0.001, retry_hhblits_with_memory_limit_if_fail=False, hhr_file=None):
    a3m_file = pathlib.Path("/tmp/")/(next(tempfile._get_candidate_names())+'.a3m')
    get_a3m_from_hhblits(input_file, a3m_file, database, maxfilt=maxfilt, realign_max=realign_max, B=B, Z=Z, n=n, e=e, retry_hhblits_with_memory_limit_if_fail=retry_hhblits_with_memory_limit_if_fail, hhr_file=hhr_file)
    call_reformat(a3m_file, output_file)


def sequence_to_msa(sequence_file, database, output_msa_file):
    get_fasta_from_hhblits(sequence_file, output_msa_file, database)


def sequence_to_processed_msa(sequence_file, database, output_msa_file, seq_id_threshold=0.8, max_nb_sequences=None):
    sequence_to_msa(sequence_file, database, output_msa_file)
    msa_processing.process_msa_for_inference(output_msa_file, output_msa_file, seq_id_threshold=seq_id_threshold, max_nb_sequences=max_nb_sequences)
