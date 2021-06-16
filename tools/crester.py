#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Description: Prepare/process structures to/from CREST
# Last update: 16-06-2021

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
  '-movable' and '-constr'. The constraint force is tunned with '-force'.

  If a bond is broken because of one side is left out of the calculation, the
  missing atom will be replaced with an "H" atom at 1Å that will be kept constrained.

  When a result structure from CREST is provided with '-x', it is inserted in its
  original position in the reference structure. The constrained atoms are used to
  align the fragment back.

"""


import os
import sys
import argparse
from copy import deepcopy
from textwrap import dedent

import numpy as np
from pdb4all import PDB, Ptable
from jacques.dynnconfig import DynnConfig
from parmed import Atom
from parmed.structure import Structure


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
    structure.translate_residues('opls')
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
    parser.add_argument('-x', metavar='.xyz', type=str, required=False,
                        help='coordinates file resulting from CREST, controls mode to/from')
    parser.add_argument('-o', metavar='<name>', type=str, required=False,
                        help='basename for output files, default taken from .dynn\n'+
                             '  w/o .xyz: <name>.xyz & <name>.constr\n'+
                             '  w/  .xyz: <name>.pdb / <name>.crd')
    parser.add_argument('-movable', metavar='<>', type=str, required=False, default='MOV',
                        help='name of atom selection group to keep movable to generate conformers (def: MOV)')
    parser.add_argument('-constr', metavar='<>', type=str, required=False, default='CONSTR',
                        help='name of atom selection group to include but keep constrained (def: CONSTR)')
    parser.add_argument('-force', metavar='#', type=float, required=False, default=0.5,
                        help='harmonic force to keep constrined atoms [Hartree*Bohr^-2] (def: 0.5)')
    args = parser.parse_args()

    ref_file     = args.c
    dynn_file    = args.f
    dynn_name    = os.path.splitext(dynn_file)[0]
    xyz_file     = args.x
    basename     = args.o or dynn_name
    movable_name = args.movable
    constr_name  = args.constr
    force        = args.force

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

    print("## CRESTER")
    if not xyz_file:
        print("# Preparing input files for CREST\n")

        # add missing atoms as H to get close shell when broken boundary bonds
        dist = 1.0  # distance to put the new H atom
        ref_pmd = pdb4all2parmed(ref_obj)
        ref_pmd_all = pdb4all2parmed(ref_obj_all)
        for n, n_all in zip(_atom_ids(ref_obj, constr_set), _atom_ids(ref_obj_all, constr_set)):
            bonds_ori = ref_pmd.atoms[n].bond_partners
            bonds_all = ref_pmd_all.atoms[n_all].bond_partners
            if len(bonds_ori) > len(bonds_all):
                idx_extra = {i.idx for i in bonds_ori} - {i.idx for i in bonds_all}
                for idx in idx_extra:
                    recept_a = ref_obj.pdb[n]
                    bonded_a = ref_obj.pdb[idx]
                    coord_recept = ref_pmd.coordinates[n]
                    coord_bonded = ref_pmd.coordinates[idx]
                    vector_bond = (coord_bonded-coord_recept) / np.linalg.norm(coord_bonded-coord_recept)
                    new_a = PDB.atom_empty
                    new_a.update({'element':"H", 'name':"H", 'segment':"X", 'resSeq':bonded_a['resSeq']})
                    new_a['x'], new_a['y'], new_a['z'] = vector_bond * dist + coord_recept
                    ref_obj_all.pdb.append(new_a)
                    constr_set.add((new_a['segment'], new_a['resSeq'], new_a['name']))
                    print(f"Added atom: 'H' -> '{bonded_a['name']}' @ "+
                          f"//{recept_a['segment']}//{recept_a['resName']}`{recept_a['resSeq']}/{recept_a['name']}")

        print("\n             movable     constr      total")
        print(f"No atoms: {len(movable_set):>10d} {len(constr_set):>10d} {ref_obj_all.natoms:>10d}\n")

        xyz_file = basename+".xyz"

        # write constrained atoms file
        constr_str = f"""\
                         $constrain
                           atoms: {",".join(map(str, _atom_ids(ref_obj_all, constr_set, shift=1)))}
                           force constant={force}
                           reference={xyz_file}
                         $end
                      """

        # write output files
        print(f"Writing constraint parameters -> '{basename}.constr'")
        with open(basename+".constr", "w") as f:
            f.write(dedent(constr_str))
        print(f"Writing coordinates -> '{basename}.xyz'")
        ref_obj_all.write_xyz(xyz_file)

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
