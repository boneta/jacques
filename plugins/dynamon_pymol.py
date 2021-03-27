#!/usr/bin/env python3

"""
  PyMOL plug-in with DYNAMON functionalities
  ==========================================

  Add extended capabilities to the PyMOL molecular viewer
  to manage files in fDynamo/DYNAMON formats

  - Write QM atoms and NOFIX residues from selection:
        write_qm  filename [, selection ]
        write_nofix  filename [, selection ]

  - Load/read extended file formats:  load  filename [...]
        .crd   -  fDynamo's coordinates file
        .dynn  -  DYNAMON options/selection file


  by Sergio Boneta Martinez
  GPLv3 2021

"""

__version__ = '0.1'

# PDB Strict formatting
# ATOM/HETATM  serial  name   altLoc  resName  chainID  resSeq  iCode  x       y     z      occupancy  tempFactor  segment  element  charge
# 1-6          7-11    13-16  17      18-20    22       23-26   27     31-38  39-46  47-54  55-60      61-66       73-76    77-78    79-80


##  DEPENDENCIES  #####################################################

import os
import fileinput

from pymol import cmd, importing
from chempy import atomic_number


##  FUNCTIONS  ########################################################

def write_sele(section, selection, filename, resolution='atom'):
    """
        Append selection to file and overwrite section if existing

        Parameters
        ----------
        section : {'QM', 'NOFIX'}
            name of section
        selection : str
            name of a PyMOL selection object
        filename : str
            file to create/append and write selection
        resolution : {'atom', 'residue', 'subsystem'}, optional
            minimum entity size to treat not whole at writting (def: 'atom')
    """

    #TODO: make it clever

    section = section.upper()

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

    # read file and delete section if existing
    try:
        in_section = False
        for line in fileinput.input(filename, inplace=True):
            if line.strip() == section:
                in_section = not in_section
            elif not in_section:
                print(line, end='')
    except FileNotFoundError:
        pass

    # write new section
    with open(filename, 'at+') as f:
        f.write(f"\n{section}\n")
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
        f.write(f"{section}\n")

    print(f" DYNAMON: {section} written to \"{os.path.abspath(filename)}\"")


def load_dynn(filename):
    """
        Load a DYNAMON selection file and read QM/NOFIX to sele objects

        Parameters
        ----------
        filename : str
            file path
    """

    print(f" DYNAMON: reading \"{filename}\"")

    # read file as list of strings
    with open(filename, "rt") as f:
        dynn_file = f.readlines()
        # remove empty lines and space-split
        dynn_file = [line.split() for line in dynn_file if line.strip()]

    i = 0
    while i < len(dynn_file):
        if dynn_file[i][0].upper() in ("QM", "NOFIX"):
            section = dynn_file[i][0].upper()
            i += 1
            # read whole selection to dict of dict of list
            sele = dict()
            segi = ""
            resi = ""
            while dynn_file[i][0].upper() != section and i < len(dynn_file):
                subsect = dynn_file[i][0].upper()
                select  = dynn_file[i][1].upper()
                if subsect == "S":
                    segi = select
                    resi = ""
                    sele[select] = dict()
                elif subsect == "R":
                    resi = select
                    sele[segi][select] = []
                elif subsect == "A":
                    sele[segi][resi].append(select)
                i += 1

            # build selection algebra
            sele1_list = []
            for segi, resdict in sele.items():
                sele1 = f"segi {segi}"
                if resdict:
                    sele1 += " & ("
                    sele2_list = []
                    for resi, namelist in resdict.items():
                        sele2 = f"resi {resi}"
                        if namelist:
                            sele2 += " & name "+"+".join(namelist)
                        sele2_list.append(sele2)
                    sele1 += " | ".join(sele2_list)+" )"
                sele1_list.append(sele1)
            sele_final = " | ".join(sele1_list)

            # selection command
            cmd.select(section, sele_final, 0, 1)
            natoms = cmd.count_atoms(section)
            print(f" DYNAMON: selection \"{section}\" defined with {natoms} atoms.")
        i += 1


def load_crd(filename, object=''):
    """
        Load a fDynamo coordinates file (.crd)

        Parameters
        ----------
        filename : str
            file path
        object : str, optional
            name of the object (def: filename prefix)
    """

    if not object:
        object = "".join(os.path.basename(filename).rpartition('.')[:-2])

    # read file as list of strings
    with open(filename, "rt") as f:
        crd_file = f.readlines()

    # convert atoms to pdb format
    atomic_number_inv = { n:elem for elem, n in atomic_number.items() }
    pdb_file = []
    a = { 'ATOM'       : "ATOM",
          'serial'     : 0,
          'name'       : "",
          'altLoc'     : "",
          'resName'    : "",
          'chainID'    : "",
          'resSeq'     : 0,
          'iCode'      : "",
          'x'          : 0.0,
          'y'          : 0.0,
          'z'          : 0.0,
          'occupancy'  : 0.0,
          'tempFactor' : 0.0,
          'segment'    : "",
          'element'    : "",
          'charge'     : "" }

    for line in crd_file:
        if line.startswith("!"): continue
        line = line.split("!")[0].split()   # to list and remove trailing comments

        if line[0].lower() == "subsystem":
            a['segment'] = str(line[2])
        elif line[0].lower() == "residue":
            a['resSeq']  = str(line[1])
            a['resName'] = str(line[2])
        elif len(line) != 6:
            continue
        else:
            a['serial']  = int(line[0])
            a['name']    = str(line[1])
            a['element'] = atomic_number_inv[int(line[2])]
            a['x']       = float(line[3])
            a['y']       = float(line[4])
            a['z']       = float(line[5])
            if a['name'] == "OXT": continue     # ignore OXT atom with usually weird coordinates
            if len(a['name']) == 3: a['name'] = ' '+a['name']     # correct alignment of atom name
            # format pdb
            formatted_line = "{:<6s}{:>5d} {:^4s}{:1s}{:>3s} {:1s}{:>4.4}{:1s}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}      {:<4s}{:>2s}{:<2s}\n" \
                .format(a['ATOM'], a['serial'], a['name'], a['altLoc'], a['resName'], a['chainID'], str(a['resSeq']), a['iCode'],
                        a['x'], a['y'], a['z'], a['occupancy'], a['tempFactor'], a['segment'], a['element'], a['charge'])
            pdb_file.append(formatted_line)

    # load as pdb string
    pdb_whole = "".join(pdb_file)
    cmd.read_pdbstr(pdb_whole, object)

    print(f" DYNAMON: \"{filename}\" loaded as \"{object}\"")


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
            .*
                default supported by PyMOL
        """

        if not format:
            format = os.path.basename(filename).rpartition('.')[-1]

        # DYNAMON dynn options/selection
        if format.lower() == "dynn":
            load_dynn(filename)

        # fDynamo's crd coordinates
        elif format.lower() == "crd":
            load_crd(filename, object)

        # original load function
        else:
            importing.load(filename, object, state, format, finish,
                           discrete, quiet, multiplex, zoom, partial)


##  WRAPPERS  #########################################################

def write_qm(dynn_file, selection="sele"):
    """Write QM atom selection to file"""
    write_sele("QM", selection, dynn_file, resolution='atom')

def write_nofix(dynn_file, selection="sele"):
    """Write NOFIX residue selection to file"""
    write_sele("NOFIX", selection, dynn_file, resolution='residue')


##  PYMOL FUNCTIONS  ##################################################
# add functions to PyMOL run environment
cmd.extend("write_qm", write_qm)
cmd.extend("write_nofix", write_nofix)
cmd.load = load_ext
cmd.extend("load", cmd.load)
