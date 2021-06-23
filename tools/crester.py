#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Description: Prepare/process structures to/from CREST
# Last update: 23-06-2021

# CREST: Conformer-Rotamer Ensemble Sampling Tool based on the xtb Semiempirical Extended Tight-Binding Program Package
# https://github.com/grimme-lab/crest / https://xtb-docs.readthedocs.io/en/latest/crest.html


"""
  Prepare/process structures to/from CREST
  ----------------------------------------

  Tested with CREST v2.11 / xTB v6.4.1
  xtb-docs.readthedocs.io/en/latest/crest.html

  Explore the conformer-rotamer lanscape of a sub-section of a bigger structure.
  This reference structure must be provided with '-c' and be in .crd/.pdb format.

  The selection of atoms is taken from a DYNAMON format compliant file always
  provided after '-f'. In this file, two selections at atomic resolution must be
  present:
    - movable atoms which conformers will be sampled
    - constrained atoms in an harmonic potential

  In case of an atom being in both groups, the atom will be included in the
  movable selection. The selection can be custom named and proper indicated with
  '-movable' and '-constr'. The constraint force is tunned with '-fc'.

  If a bond is broken because of one side is left out of any selection and it was
  part of a standard residue, the missing atom will be replaced with an "H" atom
  at 1Å that will be kept constrained.

  Several files are prepared:
    - .xyz : coordinates selected and also serve of reference for the constraints
    - .constr : constraints in xTB selection (for input after '--cinp')
    - .tmol : input for CREST in TURBOMOLE format

  Electrostatic embedding can be set by providing a force field for the charges
  of every atom and a selection of atoms after '-electr' option. 'ALL' is
  reserved as special name to include all the atoms' charges. The point charges
  hardness will be treated accordingly to the element but can be fixed with a
  larger value with '-hardness'. The corresponding embedding command is added to
  the .tmol file and the necessary .pc with the charges will is created.

  When a output structure from CREST is provided with '-x', the program changes
  its behaviour to process the results and insert it in its original position
  in the reference .pdb/.crd structure. The constrained atoms are used to align
  the fragment back.

"""


import os
import sys
import argparse
from copy import deepcopy
from textwrap import dedent

import numpy as np
from pdb4all import PDB, Ptable, aa
from jacques.dynnconfig import DynnConfig
from jacques.dyntopol import DynTopol
from parmed import Atom
from parmed.structure import Structure

# Bohr radius / Ångstrom
BOHR2A = 0.5291772109
A2BOHR = 1./BOHR2A

def _sele2set(selection: dict) -> set:
    sele_set = set()
    for segi, resis in selection.items():
        for resi, atoms in resis.items():
            for name in atoms:
                sele_set.add((segi, resi, name))
    return sele_set

def _atom_ids(structure: 'PDB', set_list: set, shift:int=0) -> list:
    return [n+shift for n, a in enumerate(structure.pdb)
            if (a['segment'], a['resSeq'], a['name']) in set_list]

def _filter_structure(structure: 'PDB', atom_ids: list) -> 'PDB':
    structure_new = deepcopy(structure)
    structure_new.pdb = [structure.pdb[n] for n in atom_ids]
    return structure_new

def pdb4all2parmed(structure:'PDB') -> 'Structure':
    structure.guess_elements()
    struc_pmd = Structure()
    for a in structure.pdb:
        atom = Atom(
                    atomic_number=Ptable[a['element']]['N'],
                    name=a['name'],
                    number=a['serial']
                    )
        struc_pmd.add_atom(
                           atom=atom,
                           resname=a['resName'],
                           resnum=a['resSeq'],
                           chain=a['chainID'],
                           segid=a['segment']
                           )
    struc_pmd.coordinates = structure.xyz
    struc_pmd.assign_bonds()
    return struc_pmd

def int_list_condensate(num_list:list, connector='-') -> list:
    """[1, 2, 2, 5, 7, 8, 9] -> ['1-2', '5', '7-9']"""
    num_list = list(sorted(set(num_list)))
    tmp = [num_list[0]]
    tmp_list = []
    for num in num_list[1:]:
        if num - 1 in tmp:
            tmp.append(num)
        else:
            tmp_list.append(deepcopy(tmp))
            tmp = [num]
    tmp_list.append(tmp)
    return [ str(nums[0]) if len(nums) == 1 else f"{nums[0]}{connector}{nums[-1]}" for nums in tmp_list]

class KabschFit:
    """
        Kabsch algorithm to fit structures minimizing RMSD

        Based on github.com/charnley/rmsd
    """

    def __init__(self, mov:np.ndarray, ref:np.ndarray):
        # centroids
        self.centroid_ref = np.mean(ref, axis=0)
        self.centroid_mov = np.mean(mov, axis=0)
        ref -= self.centroid_ref
        mov -= self.centroid_mov
        # rotation matrix based on covariance
        covar_m = np.dot(mov.T, ref)
        v, s, w = np.linalg.svd(covar_m)
        if (np.linalg.det(v) * np.linalg.det(w)) < 0.0:
            s[-1] = -s[-1]
            v[:, -1] = -v[:, -1]
        self.rot_matrix = np.dot(v, w)
        # RMSD
        diff = np.dot(mov, self.rot_matrix) - ref
        self.rmsd = np.sqrt((diff*diff).sum() / len(mov))

    def transform(self, m:np.ndarray) -> np.ndarray:
        return np.dot(m - self.centroid_mov, self.rot_matrix) + self.centroid_ref


if __name__ == '__main__':

    # parser
    parser = argparse.ArgumentParser(prog='crester.py', description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-c', metavar='.crd/.pdb', type=str, required=True,
                        help='reference coordinates file')
    parser.add_argument('-f', metavar='.dynn', type=str, required=True,
                        help='configuration file in DYNAMON format to take atomic selections and constraints')
    parser.add_argument('-ff', metavar='.ff', type=str, required=False,
                        help='force field in fDynamo format to take electrostatic charges')
    parser.add_argument('-x', metavar='.xyz', type=str, required=False,
                        help='coordinates file resulting from CREST, controls mode to/from')
    parser.add_argument('-o', metavar='<name>', type=str, required=False,
                        help='basename for output files, default taken from .dynn\n'+
                             '  w/o .xyz: .xyz, .tmol, .constr [, .pc]\n'+
                             '  w/  .xyz: .pdb/.crd')
    parser.add_argument('-movable', metavar='<>', type=str, required=False, default='MOV',
                        help='name of atom selection group to keep movable to generate conformers (def: MOV)')
    parser.add_argument('-constr', metavar='<>', type=str, required=False, default='CONSTR',
                        help='name of atom selection group to include but keep constrained (def: CONSTR)')
    parser.add_argument('-electr', metavar='<>', type=str, required=False,
                        help='name of atom selection group to include as point charges (special name: ALL)')
    parser.add_argument('-fc', metavar='#', type=float, required=False, default=0.1,
                        help='harmonic force to keep constrined atoms [Hartree*Bohr^-2] (def: 0.1)')
    parser.add_argument('-hardness', metavar='#', type=int, required=False,
                        help='chemical hardness behaviour for partial charges (def: elemental)')
    args = parser.parse_args()

    ref_file     = args.c
    dynn_file    = args.f
    dynn_name    = os.path.splitext(dynn_file)[0]
    ff_file      = args.ff
    xyz_file     = args.x
    basename     = args.o or dynn_name
    movable_name = args.movable
    constr_name  = args.constr
    electr_name  = args.electr
    fc           = args.fc
    hardness     = args.hardness
    electr_flg   = ff_file and electr_name


    # initalizate objects & read files
    ref_obj = PDB()
    ref_extension = os.path.splitext(ref_file)[-1].lower()
    if ref_extension == '.pdb':
        ref_obj.read(ref_file)
    elif ref_extension == '.crd':
        ref_obj.read_crd(ref_file)
    else:
        sys.exit(f"ERROR: Unknown reference file extension: '{ref_extension}'")
    dynn_obj = DynnConfig(dynn_file)

    #TODO: add special constraints from dynn

    # sets of tuples of (segment, residue, atom)
    movable_set = _sele2set(dynn_obj.selection.get(movable_name, dict()))
    constr_set  = _sele2set(dynn_obj.selection.get(constr_name, dict())) - movable_set
    all_set     = movable_set | constr_set

    # check selections sizes
    if not (movable_set and constr_set):
        sys.exit("ERROR: Missing a selection atoms")
    elif len(constr_set) < 3:
        sys.exit("ERROR: 3+ constrained atoms are required")

    # create objects with only atoms based on selections
    ref_obj_all    = _filter_structure(ref_obj, _atom_ids(ref_obj, all_set))
    ref_obj_constr = _filter_structure(ref_obj, _atom_ids(ref_obj, constr_set))

    if not xyz_file:
        print("# Preparing input files for CREST\n")

        # add missing atoms as H to get close shell when broken boundary bonds
        new_dist = 1.0  # distance to put the new H atom
        ref_pmd = pdb4all2parmed(ref_obj)
        ref_pmd_all = pdb4all2parmed(ref_obj_all)
        coord_all = np.array(ref_obj_all.xyz)
        for n, n_all in zip(_atom_ids(ref_obj, constr_set), _atom_ids(ref_obj_all, constr_set)):
            bonds_ori = ref_pmd.atoms[n].bond_partners
            bonds_all = ref_pmd_all.atoms[n_all].bond_partners
            if len(bonds_ori) > len(bonds_all):
                idx_extra = {i.idx for i in bonds_ori} - {i.idx for i in bonds_all}
                for idx in idx_extra:
                    recept_a = ref_obj.pdb[n]
                    bonded_a = ref_obj.pdb[idx]
                    if recept_a['resName'] not in aa: break     # skip if not a residue
                    coord_recept = ref_pmd.coordinates[n]
                    coord_bonded = ref_pmd.coordinates[idx]
                    vector_bond = (coord_bonded-coord_recept) / np.linalg.norm(coord_bonded-coord_recept)
                    coord_new = vector_bond * new_dist + coord_recept
                    # check if new coordinates are too close to any other atom
                    diff = np.delete(coord_all, n_all, 0) - coord_new
                    dist = np.sqrt((diff*diff).sum(axis=1))
                    if np.amin(dist) < 1.1: continue
                    # add new atom
                    new_a = deepcopy(PDB.atom_empty)
                    new_a.update({'element':"H", 'name':"H", 'segment':"X", 'resSeq':bonded_a['resSeq']})
                    new_a['x'], new_a['y'], new_a['z'] = coord_new
                    ref_obj_all.pdb.append(new_a)
                    constr_set.add((new_a['segment'], new_a['resSeq'], new_a['name']))
                    print(f"Boundary atom to H: {PDB.pymol_sele_macro(bonded_a)} connected to {PDB.pymol_sele_macro(recept_a)}")

        print("\n             movable     constr      total")
        print(f"No atoms: {len(movable_set):>10d} {len(constr_set):>10d} {ref_obj_all.natoms:>10d}\n")

        # constrained atoms
        constr_str = f"""\
                      $constrain
                          atoms: {",".join(int_list_condensate(_atom_ids(ref_obj_all, constr_set, shift=1)))}
                          force constant={fc}
                          reference={basename}.xyz
                      """
        constr_str = dedent(constr_str)

        # electrostatic embedding
        if electr_flg:
            electr_set = {(a['segment'], a['resSeq'], a['name']) for a in ref_obj.pdb} if electr_name.upper() == "ALL" \
                         else _sele2set(dynn_obj.selection.get(electr_name, dict()))
            electr_set -= all_set
            if not electr_set:
                sys.exit("ERROR: Electrostatic selection empty or not found")
            ref_obj_electr = _filter_structure(ref_obj, _atom_ids(ref_obj, electr_set))
            # add charge property
            topol = DynTopol(ff_file)
            for atom in ref_obj_electr.pdb:
                try:
                    atom['charge'] = topol.top['residues'][atom['resName']]['atoms'][atom['name']]['charge']
                except KeyError:
                    print(f"WARNING: Charge not found for atom '{PDB.pymol_sele_macro(atom)}'")
                    atom['charge'] = 0.
            # change hardness if requested
            if hardness:
                ref_obj_electr.clean_field('element', str(hardness))
            # add embedding command
            constr_str += f"$embedding\n    input={basename}.pc\n"
            # write point charges file
            ref_obj_electr.remove('charge', 0.)     # remove 0 charges
            print(f"Writing point charges -> '{basename}.pc'")
            with open(basename+".pc", "w") as f:
                # write coordinates
                f.write(f"{ref_obj_electr.natoms}\n")
                for atom in ref_obj_electr.pdb:
                    f.write(" {:>18.10f}  {:>18.10f} {:>18.10f} {:>18.10f}       {:<6s}\n"
                            .format(atom['charge'], atom['x']*A2BOHR, atom['y']*A2BOHR, atom['z']*A2BOHR, atom['element']))

        # XYZ
        print(f"Writing coordinates -> '{basename}.xyz'")
        ref_obj_all.write_xyz(basename+".xyz")
        # CONSTRAINTS
        print(f"Writing constraints -> '{basename}.constr'")
        with open(basename+".constr", "w") as f:
            f.write(dedent(constr_str))
            f.write("$end\n")
        # TURBOMOLE
        print(f"Writing CREST input -> '{basename}.tmol'")
        with open(basename+".tmol", "w") as f:
            # write constraints (not recognized by CREST)
            f.write(dedent(constr_str))
            # write coordinates
            f.write("$coord\n")
            for atom in ref_obj_all.pdb:
                f.write(" {:>18.10f} {:>18.10f} {:>18.10f}       {:<6s}\n"
                        .format(atom['x']*A2BOHR, atom['y']*A2BOHR, atom['z']*A2BOHR, atom['element'].lower()))
            f.write("$end\n")

    else:
        print("# Processing output from CREST")

        xyz_obj = PDB()
        xyz_obj.read_xyz(xyz_file)
        xyz_obj_constr = _filter_structure(xyz_obj, _atom_ids(ref_obj_all, constr_set))

        if ref_obj_all.natoms != xyz_obj.natoms:
            sys.exit("ERROR: Missmatching number of selected atoms between .xyz and reference structure")

        # align xyz to reference constrained coordinates
        kabsch = KabschFit(xyz_obj_constr.xyz, ref_obj_constr.xyz)
        print(f"RMSD of fitted constrained atoms [A]:  {kabsch.rmsd:.5f}")

        # transplant transformed coordinates from xyz to reference crd
        transformed_coord = kabsch.transform(xyz_obj.xyz)
        for n_crd, n_xyz in enumerate(_atom_ids(ref_obj, all_set)):
            ref_obj.pdb[n_crd]['x'], \
            ref_obj.pdb[n_crd]['y'], \
            ref_obj.pdb[n_crd]['z'] = map(float, transformed_coord[n_xyz])

        # write output files
        print(f"Writing coordinates -> '{basename}{ref_extension}'")
        if ref_extension == '.pdb':
            ref_obj.write(basename+".pdb")
        elif ref_extension == '.crd':
            ref_obj.write_crd(basename+".crd")
