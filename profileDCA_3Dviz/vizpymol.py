import argparse, sys, pathlib
import pymol

from profileDCA_utils import io_management as iom
from profileDCA_3Dviz import contacts_management, pdb_utils


def script_init(fout, pdb_id=None, pdbfile=None):
    if pdbfile is not None:
        fout.write(f'load {pdbfile}\n')
    elif pdb_id is not None:
        fout.write(f'fetch {pdbfile}\n')
    else:
        raise Exception("must provide pdb file or pdb id")
    fout.write(f'show cartoon\n')
    fout.write(f'hide lines\n')
    fout.write(f'hide nonbonded\n')
    fout.write(f'set cartoon_color, grey\n')
    fout.write(f'set cartoon_transparency, 0.2\n')

    fout.write(f'set_color contact_coupling_color1, [0.00 , 0.53 , 0.22] # green\n')
    fout.write(f'set_color distant_coupling_color1, [1.00 , 0.00 , 0.00] # red\n')

    fout.write(f'set_color contact_coupling_color2, [0.0 , 0.0 , 1.0] # blue\n')
    fout.write(f'set_color distant_coupling_color2, [1.0 , 1.0 , 0.0] # yellow\n')

    fout.write(f'set_color contact_coupling_color3, [0.0 , 0.75 , 0.75] # teal\n')
    fout.write(f'set_color distant_coupling_color3, [1.0 , 0.5 , 0.0] # orange\n')


def launch_pymol(fout, pdb_id=None, pdbfile=None):
    script_init(fout, pdb_id, pdbfile)
    pymol.cmd.reinitialize()
    # pymol.finish_launching(['pymol'])
    if pdbfile is not None:
        pymol.cmd.load(pdbfile)
    elif pdb_id is not None:
        pymol.cmd.fetch(pdb_id, async_=0, type="pdb")
    else:
        raise Exception("must provide pdb file or pdb id")
    pymol.cmd.show('cartoon')
    pymol.cmd.hide('lines')
    pymol.cmd.hide('nonbonded')
    pymol.cmd.set('cartoon_color', 'grey')
    pymol.cmd.set('cartoon_transparency', 0.2)

    pymol.cmd.set_color(
        'contact_coupling_color1', [0.00, 0.53, 0.22])
    pymol.cmd.set_color(
        'distant_coupling_color1', [1.00, 0.00, 0.00])
    pymol.cmd.set_color(
        'contact_coupling_color2', [0.0, 0.0, 1.0])
    pymol.cmd.set_color(
        'distant_coupling_color2', [1.0, 1.0, 0.0])
    pymol.cmd.set_color(
        'contact_coupling_color3', [0.0, 0.75, 0.75])
    pymol.cmd.set_color(
        'distant_coupling_color3', [1.0, 0.5, 0.0])


def script_coupling(fout, pdb_coupling, strength, color, chain_id='A'):
    pos1 = pdb_coupling[0]
    pos2 = pdb_coupling[1]
    fout.write(f'bond {chain_id}/{pos1}/CA, {chain_id}/{pos2}/CA\n')
    fout.write(f'select {chain_id}/{pos1}+{pos2}/CA\n')
    fout.write(f'set_bond stick_color, {color}, sele\n')
    fout.write(f'set_bond stick_radius, {strength}, sele\n')
    fout.write(f'show sticks, sele\n')
    fout.write(f'label sele, resi\n')


def show_coupling(fout, pdb_coupling, strength, color, chain_id='A'):
    script_coupling(fout, pdb_coupling, strength, color, chain_id)
    pos1 = pdb_coupling[0]
    pos2 = pdb_coupling[1]
    pymol.cmd.select("coupling", "resi "+str(pos1)+" and name CA and chain " +
                     chain_id+" + resi "+str(pos2)+" and name CA and chain "+chain_id)
    pymol.cmd.bond("resi "+str(pos1)+" and name CA and chain " +
                   chain_id, "resi "+str(pos2)+" and name CA and chain "+chain_id)
    pymol.cmd.set_bond("stick_color", color, "coupling")
    pymol.cmd.set_bond("stick_radius", strength, "coupling")
    pymol.cmd.show('sticks', "coupling")
    pymol.cmd.label("coupling", 'resi')


def show_n_couplings(fout, nb_couplings, pdb_seq_couplings_dict, pdb_chain, coupling_sep_min=2, thickness=1, colors={True: 'contact_coupling_color1', False: 'distant_coupling_color1'}, contact_distance=8):
    #pdb_chain = fm.get_pdb_chain(pdb_id, pdb_file, chain_id)
    n = 0
    for i, (c, score) in enumerate(pdb_seq_couplings_dict.items()):
        if n < nb_couplings:
            if abs(c[0]-c[1]) > coupling_sep_min:
                strength = score*thickness
                show_coupling(fout, c, strength,
                              colors[pdb_utils.is_pdb_pair_contact(*c, pdb_chain, contact_threshold=contact_distance)])
                n += 1


def show_predicted_contacts_with_pymol(fout, potts_folders, pdb_id=None, chain_id='A', pdb_file=None, top=20, coupling_sep_min=3, thickness=1, auto_top=False, wij_cutoff=None, normalize=False, debug_mode=False, contact_distance=8, **kwargs):

    if (pdb_file is None) and (pdb_id is not None):
        name = potts_folders[0]/pdb_id
        pdb_file = pdb_utils.fetch_pdb_file(pdb_id, name)
    pdb_chain = pdb_utils.get_pdb_chain(pdb_file, pdb_id, chain_id)

    pdb_couplings_dicts = []
    tops = []
    for potts_folder in potts_folders:
        mrf = iom.mrf_from_folder(potts_folder)
        original_sequence = iom.get_first_sequence_in_fasta_file(potts_folder/"sequence.fasta")
        pdb_couplings_dict = contacts_management.get_contact_scores_with_pdb_indexes(mrf, original_sequence, pdb_chain)
        if wij_cutoff:
            cutindex = contacts_management.get_cutoff_smaller_than(pdb_couplings_dict, wij_cutoff)
            pdb_couplings_dict = contacts_management.get_smaller_dict(pdb_couplings_dict, cutindex)
            nb_couplings = len(pdb_couplings_dict)
        else:
            nb_couplings = top
        pdb_couplings_dict = contacts_management.remove_couplings_too_close(pdb_couplings_dict, coupling_sep_min)
        pdb_couplings_dicts.append(pdb_couplings_dict)
        tops.append(nb_couplings)

    launch_pymol(fout, pdb_id, pdb_file)

    exclus_overlap = contacts_management.get_exclus_overlaps(pdb_couplings_dicts, tops)

    if normalize:
        for k in range(len(exclus_overlap)):
            exclus_overlap[k] = get_normalized_ordered_dict(exclus_overlap[k])

    for d, colors in zip(exclus_overlap, [{True: 'contact_coupling_color1', False: 'distant_coupling_color1'}, {True: 'contact_coupling_color2', False: 'distant_coupling_color2'}, {True: 'contact_coupling_color3', False: 'distant_coupling_color3'}]):
        show_n_couplings(fout, len(d), d, pdb_chain,
                         coupling_sep_min=coupling_sep_min, thickness=thickness, colors=colors, contact_distance=contact_distance)


def main(args=sys.argv[1:]):
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--potts_folders', help="Potts folder(s)",
                        type=pathlib.Path, nargs='+', required=True)
    parser.add_argument('--pdb_file', help="PDB file",
                        type=pathlib.Path, default=None)
    parser.add_argument('-id', '--pdb_id', help="PDB id")
    parser.add_argument('-cid', '--chain_id',
                        help="PDB chain id (default : A)", default='A')
    parser.add_argument('-sep', '--coupling_sep_min',
                        help="Min. nb residues between members of a coupling (default : 3)", default=3, type=int)
    parser.add_argument('--contact_distance',
            help="Distance threshold in Angstrom for a contact (default: 8)", default=8, type=float)
    parser.add_argument(
        '-n', '--top', help="Nb of couplings displayed (default : 20)", type=int, default=20)
    parser.add_argument(
        '--wij_cutoff', help="||wij|| <= wij_cutoff are removed (default : None)", default=None, type=float)
    parser.add_argument('--auto_top', help="Nb couplings displayed = elbow of the score curve (default : False)",
                        default=False, action='store_true')
    parser.add_argument(
        '-t', '--thickness', help="Couplings thickness factor (default : 1)", type=float, default=1)
    parser.add_argument('--normalize', help="Normalize coupling values (default : don't normalize)",
                        default=False, action='store_true')
    parser.add_argument('--debug_mode', help=argparse.SUPPRESS,
                        default=False, action='store_true')
    parser.add_argument('--out_session_file', '-pse',
                        help="PyMOL output session file (must end in .pse) (default : /tmp/tmp_pymol_session_file.pse)",
                        type=pathlib.Path, default=pathlib.Path('/tmp/tmp_pymol_session_file.pse'))
    parser.add_argument('--out_script_file', '-pml',
                        help="PyMOL output script file (must end in .pml) (default : /tmp/tmp_pymol_script_file.pml)",
                        type=pathlib.Path, default=pathlib.Path('/tmp/tmp_pymol_session_file.pml'))

    args = vars(parser.parse_args(args))

    with open(args['out_script_file'], 'w') as fout:
        show_predicted_contacts_with_pymol(fout, **args)
    pymol.cmd.save(str(args["out_session_file"]))
    print("PyMOL session saved at "+str(args["out_session_file"]))
    print("PyMOL script saved at "+str(args["out_script_file"]))


if __name__ == "__main__":
    main()
