#!/usr/bin/env python3

"""
  PyMOL plug-in with DYNAMON functionalities
  ==========================================

  Manage atom/residue selections in fDynamo/DYNAMON format

  - Write QM atoms and NOFIX residues from selection:
        write_qm  dynn_file [, selection ]
        write_nofix  dynn_file [, selection ]

  - Load/read extended file formats:  load  filename [...]
        .crd   -  fDynamo's coordinates file
        .dynn  -  DYNAMON options/selection file
"""

__version__ = '0.1'


##  DEPENDENCIES  #####################################################

import os
from pymol import cmd, importing


##  DEFAULT OPTIONS  ##################################################

_selection_def = 'sele'


##  FUNCTIONS  ########################################################

def write_sele(section, selection, dynn_file, resolution='atom'):
    """
        Write selection to file (append)

        Parameters
        ----------
        section : {'QM', 'NOFIX'}
            name of section
        selection : str
            name of a PyMOL selection object
        dynn_file : str
            file to create/append and write selection
        resolution : {'atom', 'residue', 'subsystem'}, optional
            minimum entity size to treat not whole at writting (def: 'atom')
    """

    #TODO: make it clever

    # get list of sele atoms with dictionary of properties
    natoms = 0
    atoms = []
    obj_list = cmd.get_object_list(selection)
    for obj in obj_list:
        model = cmd.get_model(selection+" and "+obj)
        natoms += model.nAtom
        for a in model.atom:
            name = str(a.name)
            resi = int(a.resi)
            segi = str(a.segi)
            atoms.append({ 'name' : name,
                           'resi' : resi,
                           'segi' : segi })

    print(f" DYNAMON: Number of atoms in \"{selection}\": {natoms}")

    # write to file
    f = open(dynn_file, 'at+')

    f.write(f"\n{section.upper()}\n")

    segis = sorted(list(set([a['segi'] for a in atoms])))   # unique subsystem list
    for segi in segis:
        f.write(" "*4+f"S {segi}\n")
        if resolution == 'subsystem': continue

        resis = sorted(list(set([a['resi'] for a in atoms if a['segi']==segi])))
        for resi in resis:
            f.write(" "*8+f"R {resi}\n")
            if resolution == 'residue': continue

            names = sorted(list(set([a['name'] for a in atoms if a['segi']==segi and a['resi']==resi])))
            for name in names:
                f.write(" "*12+f"A {name}\n")

    f.write(f"{section.upper()}\n\n")

    f.close()
    print(f" DYNAMON: {section} written to \"{os.path.abspath(dynn_file)}\"")


def load_ext(filename, object='', state=0, format='', finish=1,
             discrete=-1, quiet=1, multiplex=None, zoom=-1, partial=0,
             mimic=1, object_props=None, atom_props=None, *, _self=cmd):
        """
            Wrapper to load function with extended funtionality

            Formats
            -------
            .crd
                fDynamo's coordinates file
            .dynn
                DYNAMON options/selection file
        """

        if not format:
            file_dot_separated = os.path.basename(filename).rpartition('.')
            name   = "".join(file_dot_separated[:-2])
            format = file_dot_separated[-1]

        # DYNAMON dynn options/selection
        if format.lower() == "dynn":
            print(f" DYNAMON: reading \"{filename}\"")
            print(f" DYNAMON: Loading .dynn not yet implemented")
            return

        # fDynamo's crd coordinates
        if format.lower() == "crd":
            print(f" DYNAMON: \"{filename}\" loaded as \"{name}\"")
            print(f" DYNAMON: Loading .crd not yet implemented")
            return

        # original load function
        else:
            importing.load(filename, object, state, format, finish,
                           discrete, quiet, multiplex, zoom, partial)


##  WRAPPERS  #########################################################

def write_qm(dynn_file, selection=_selection_def):
    """Write QM atom selection to file"""
    write_sele("QM", selection, dynn_file, resolution='atom')

def write_nofix(dynn_file, selection=_selection_def):
    """Write NOFIX residue selection to file"""
    write_sele("NOFIX", selection, dynn_file, resolution='residue')


##  PYMOL FUNCTIONS  ##################################################
# add functions to PyMOL run environment
cmd.extend("write_qm", write_qm)
cmd.extend("write_nofix", write_nofix)
cmd.load = load_ext
cmd.extend("load", cmd.load)
