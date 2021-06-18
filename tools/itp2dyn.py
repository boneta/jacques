#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Description: Convert parameters from GROMACS' .itp (AMBER/GAFF) to DYNAMO's .ff (OPLS)
# Last update: 18-06-2021


import sys
import re
import readline
import argparse
import fileinput
from typing import Union

from jacques.dyntopol import DynTopol
from amber2opls import DihedralParam


##  INFO  #############################################################

"""
    From ParmEd (gromacstop.py) https://github.com/ParmEd/ParmEd

    Gromacs uses "funct" flags in its parameter files to indicate what kind of
    functional form is used for each of its different parameter types. This is
    taken from the topdirs.c source code file along with a table in the Gromacs
    user manual.

    Bonds
    -----
    1  - F_BONDS      : simple harmonic potential
    2  - F_G96BONDS   : fourth-power potential
    3  - F_MORSE      : morse potential
    4  - F_CUBICBONDS : cubic potential
    5  - F_CONNBONDS  : not even implemented in GROMACS
    6  - F_HARMONIC   : seems to be the same as (1) ??
    7  - F_FENEBONDS  : finietely-extensible-nonlinear-elastic (FENE) potential
    8  - F_TABBONDS   : bond function from tabulated function
    9  - F_TABBONDSNC : bond function from tabulated function (no exclusions)
    10 - F_RESTRBONDS : restraint bonds

    Angles
    ------
    1 - F_ANGLES            : simple harmonic potential
    2 - F_G96ANGLES         : cosine-based angle potential
    3 - F_CROSS_BOND_BONDS  : bond-bond cross term potential
    4 - F_CROSS_BOND_ANGLES : bond-angle cross term potential
    5 - F_UREY_BRADLEY      : Urey-Bradley angle-bond potential
    6 - F_QUARTIC_ANGLES    : 4th-order polynomial potential
    7 - F_TABANGLES         : angle function from tabulated function
    8 - F_LINEAR_ANGLES     : angle function from tabulated function
    9 - F_RESTRANGLES       : restricted bending potential

    Dihedrals
    ---------
    1  - F_PDIHS     : periodic proper torsion potential [ k(1+cos(n*phi-phase)) ]
    2  - F_IDIHS     : harmonic improper torsion potential
    3  - F_RBDIHS    : Ryckaert-Bellemans torsion potential
    4  - F_PIDIHS    : periodic harmonic improper torsion potential (same as 1)
    5  - F_FOURDIHS  : Fourier dihedral torsion potential
    8  - F_TABDIHS   : dihedral potential from tabulated function
    9  - F_PDIHS     : Same as 1, but can be multi-term
    10 - F_RESTRDIHS : Restricted torsion potential
    11 - F_CBTDIHS   : combined bending-torsion potential
"""


##  CONSTANTS  ########################################################

cal2J = 4.184      # cal -> J
ptable_simple = {'D':1,'T':1,'H':1,'C':6,'N':7,'O':8,'F':9,'P':15,'S':16}


##  GROMACS TOPOLOGY CLASS  ###########################################

class GMXTopology():

    #TODO: support more parameters types

    def __init__(self, file_inp=None):
        self.atomtypes = dict()
        self.atomtypes_rosetta = dict()
        self.moleculetype = dict()
        self.system = ""
        self.molecules = dict()
        self.opls = DynTopol()
        self.missing = dict()
        self.found = DynTopol()
        if file_inp is not None:
            self.file_inp = file_inp
            self.read_itp(file_inp)
        else:
            self.file_inp = ""

    @property
    def molnames(self) -> list:
        return list(self.moleculetype.keys())

    @property
    def nmolec(self) -> int:
        return len(self.moleculetype)

    def read_itp(self, file_inp:str, atomtypes:str='literal') -> None:

        # TODO: implements more parameter types

        # read topology file
        with open(file_inp, 'rt') as f:
            ftop = f.readlines()
            ftop = list(map(str.strip, ftop))  # remove trailing spaces

        def _check_num_param(param_list, section, nline, nmin, nmax):
            """Check correct number of parameters or throw an error"""
            nparam = len(param_list)
            msg = f"{section} -> Wrong nuber of parameters (file {file_inp}, line {nline+1})"
            if not nmin <= nparam <= nmax: raise ValueError(msg)

        def _read_atomtype(attype:str, mode:str='literal'):
                if mode not in ('literal','inter','element'):
                    raise ValueError("Unkown atomtype reading mode")
                if mode == 'literal':
                    readed_attype = attype
                elif mode == 'inter':
                    readline.set_startup_hook(lambda: readline.insert_text(attype.upper()))
                    try:
                        readed_attype = input("Interactive atomtype reading: ")
                    finally:
                        readline.set_startup_hook()
                elif mode == 'element':
                    readed_attype = re.findall(r'[a-zA-Z]', attype)[0]
                self.atomtypes_rosetta[attype] = readed_attype
                return readed_attype

        # process line by line
        current_sele = None
        for nline, line in enumerate(ftop):
            if not line.strip() or line.startswith((";","#")): continue

            if line.startswith('['):
                current_sele = line.split()[1]

            elif current_sele == 'atomtypes': # -----------------------
                """
                    Fields
                    ------
                    0 : attype      : str,   mandatory  (nobond_type)
                    1 : bond_type   : str,   optional
                    2 : atnum       : int,   optional   (atomic number)
                    3 : mass        : float, mandatory
                    4 : charge      : float, mandatory
                    5 : ptype       : str,   mandatory  (particle type, one letter)
                    6 : sigma       : float, mandatory
                    7 : epsilon     : float, mandatory
                """
                param = { 'bond_type' : "",
                          'atnum' : "",
                          'mass' : "",
                          'charge' : "",
                          'ptype' : "",
                          'sigma' : "",
                          'epsilon' : "" }
                words = line.split(";")[0].split()  # remove comments and list
                _check_num_param(words, current_sele, nline, 6, 8)
                if len(words) == 8:  # both optionals are present
                    param['bond_type'] = str(words[1])
                    param['atnum'] = int(words[2])
                elif len(words) == 7:  # one optional is present
                    try:
                        param['atnum'] = int(words[1])
                    except ValueError:
                        param['bond_type'] = str(words[1])
                # assign rest, from last to first
                words = words[::-1]
                attype = _read_atomtype(str(words[-1]), atomtypes)
                param['epsilon'] = float(words[0])
                param['sigma'] = float(words[1])
                param['ptype'] = str(words[2])
                param['charge'] = float(words[3])
                param['mass'] = float(words[4])
                # guess atomic number if absent
                if not param['atnum']: param['atnum'] = ptable_simple[re.findall(r'[a-zA-Z]', words[-1])[0].upper()]
                self.atomtypes[attype] = param

            elif current_sele == 'moleculetype': # --------------------
                param = { 'nrexcl' : "",
                          'atoms' : dict(),
                          'bonds' : dict(),
                          'pairs' : dict(),
                          'angles' : dict(),
                          'dihedrals' : dict() }

                words = line.split(";")[0].split()  # remove comments and list
                _check_num_param(words, current_sele, nline, 2, 2)
                molname = str(words[0])
                param['nrexcl'] = int(words[1])
                if molname in self.molnames:
                    raise ValueError(f"moleculetype -> Duplicate definition of molecule ({molname}, file {file_inp}, line {nline+1})")
                self.moleculetype[molname] = param
                current_molec = molname

            elif current_sele == 'atoms': # ---------------------------
                param = { 'attype' : "",
                          'resid' : "",
                          'resname' : "",
                          'name' : "",
                          'cgnr' : "",
                          'charge' : "",
                          'mass' : "" }
                words = line.split(";")[0].split()  # remove comments and list
                _check_num_param(words, current_sele, nline, 8, 8)
                nr = int(words[0])
                param['attype'] = self.atomtypes_rosetta[str(words[1])]
                param['resid'] = int(words[2])
                param['resname'] = str(words[3])
                param['name'] = str(words[4])
                param['cgnr'] = int(words[5])
                param['charge'] = float(words[6])
                param['mass'] = float(words[7])
                self.moleculetype[current_molec]['atoms'][nr] = param

            elif current_sele == 'bonds': # ---------------------------
                param = { 'funct' : "",
                          'r' : "",
                          'k' : "" }
                words = line.split(";")[0].split()  # remove comments and list
                _check_num_param(words, current_sele, nline, 3, 5)
                atoms = tuple(int(w) for w in words[:2])
                param['funct'] = int(words[2])
                if param['funct'] != 1: raise ValueError(f"bonds -> only funct==1 allowed (file {file_inp}, line {nline+1})")
                if len(words) == 5:
                    param['r'] = float(words[3])
                    param['k'] = float(words[4])
                self.moleculetype[current_molec]['bonds'][atoms] = param

            elif current_sele == 'pairs': # ---------------------------
                param = { 'funct' : "" }
                words = line.split(";")[0].split()  # remove comments and list
                _check_num_param(words, current_sele, nline, 3, 3)
                atoms = tuple(int(w) for w in words[:2])
                param['funct'] = int(words[2])
                if param['funct'] != 1: raise ValueError(f"pairs -> only funct==1 allowed (line {nline+1})")
                self.moleculetype[current_molec]['pairs'][atoms] = param

            elif current_sele == 'angles': # --------------------------
                param = { 'funct' : "",
                          'theta' : "",
                          'cth' : "" }
                words = line.split(";")[0].split()  # remove comments and list
                _check_num_param(words, current_sele, nline, 4, 6)
                atoms = tuple(int(w) for w in words[:3])
                param['funct'] = int(words[3])
                if param['funct'] != 1: raise ValueError(f"angles -> only funct==1 allowed (file {file_inp}, line {nline+1})")
                if len(words) == 6:
                    param['theta'] = float(words[4])
                    param['cth'] = float(words[5])
                self.moleculetype[current_molec]['angles'][atoms] = param

            elif current_sele == 'dihedrals': # -----------------------
                param = { 'funct' : 0,
                          'phase' : [],
                          'kd' : [],
                          'pn' : [],
                          'opls' : [] }
                words = line.split(";")[0].split()  # remove comments and list
                _check_num_param(words, current_sele, nline, 5, 8)
                atoms = tuple(int(w) for w in words[:4])
                funct = int(words[4])
                if funct not in (1, 4, 9):
                    raise ValueError(f"dihedrals -> funct not allowed (file {file_inp}, line {nline+1})")
                elif len(words)>=8:
                    # check if dihedral already stored for this atoms
                    if atoms in self.moleculetype[current_molec]['dihedrals']:
                        param = self.moleculetype[current_molec]['dihedrals'][atoms]
                    param['phase'].append(float(words[5]))
                    param['kd'].append(float(words[6]))
                    param['pn'].append(int(words[7]))
                else:
                    param['phase'].append("")
                    param['kd'].append("")
                    param['pn'].append("")
                param['funct'] = funct
                # append to general list only if new dihedral parameter
                if len(param['phase']) == 1: self.moleculetype[current_molec]['dihedrals'][atoms] = param

            elif current_sele == 'system':  # -------------------------
                self.system = str(line.split(";")[0])  # remove comments

            elif current_sele == 'molecules':  # ----------------------
                words = line.split(";")[0].split()  # remove comments and list
                self.molecules[str(words[0])] = int(words[1])

            elif current_sele in {'cmap', 'defaults', 'settles', 'exclusions',
                                  'nonbond_params', 'bondtypes', 'angletypes',
                                  'dihedraltypes', 'cmaptypes', 'pairtypes'}:
                continue
            else:
                raise ValueError(f"Unclassifiable selection: '{current_sele}' (file {file_inp}, line {nline+1})")

    def gen_opls(self) -> None:

        self.opls.__init__()
        top = self.opls.top

        # atom types
        for attype, atom in self.atomtypes.items():
            top['atomtypes'][attype.upper()] = {'atnum':atom['atnum'], 'sigma':atom['sigma']*10,
                                                'epsilon':atom['epsilon'], 'comment':""}

        # atoms per residue and unified bonds/angles/dihedrals/impropers
        param = DihedralParam()
        for molname, molec in self.moleculetype.items():
            # atoms in residue
            top['residues'][molname] = DynTopol.empty_param('residues')
            for nr, atom in molec['atoms'].items():
                if atom['name'] in top['residues'][molname]['atoms']:
                    raise ValueError(f"opls conversion -> duplicated atom definition (residue {molname}, atom {atom['name']})")
                top['residues'][molname]['atoms'][atom['name']] = {'atomtype':atom['attype'].upper(), 'charge':atom['charge']}
            # bonds
            for atoms, bond in molec['bonds'].items():
                attypes = tuple(molec['atoms'][nr]['attype'].upper() for nr in atoms)
                top['bonds'][attypes] = {'k':bond['k']/1000, 'r':bond['r']*10, 'comment':""}
            # angles
            for atoms, angle in molec['angles'].items():
                attypes = tuple(molec['atoms'][nr]['attype'].upper() for nr in atoms)
                top['angles'][attypes] = {'k':angle['cth']/10, 'theta':angle['theta'], 'comment':""}
            # dihedrals/impropers
            for atoms, dihedral in molec['dihedrals'].items():
                divider = 2*cal2J   # fitted by comparison
                attypes = tuple(molec['atoms'][nr]['attype'].upper() for nr in atoms)
                amberlist = []
                for i in range(len(dihedral['phase'])):
                    amberlist.extend([divider, dihedral['kd'][i], dihedral['phase'][i], dihedral['pn'][i]])
                param.list2amber(amberlist, phase_units='degree')
                if param.amber2opls():
                    # dihedral : 9 / improper : 4
                    param_type = 'dihedrals' if dihedral['funct'] == 9 else 'impropers'
                    top[param_type][attypes] = {'v':param.p['opls'], 'comment':""}

        self.opls.recount_top_n()
        self.opls.remove_duplicates()

    def read_miss(self, missing:list) -> None:
        self.missing = { 'bonds' : set(),
                         'angles' : set(),
                         'dihedrals' : set(),
                         'impropers' : set() }

        # keep only error records with "Missing"
        missing = [line for line in missing if line.strip().startswith('Missing')]

        for line in missing:
            ptype = line.split()[1]
            line = line.split(":")[1].split()
            if ptype in ('bond', 'angle', 'dihedral', 'improper'):
                self.missing[ptype+'s'].add(tuple(line))
                # remove inverse if not palindromic
                if line != line[::-1]: self.missing[ptype+'s'].discard(tuple(line[::-1]))
            else:
                raise NameError(f'Unknown missing parameter type: {ptype}')

        for ptype in self.missing.keys():
            self.missing[ptype] = sorted(list(self.missing[ptype]))

    def find_miss(self) -> None:
        self.found = DynTopol()
        if not self.opls: self.gen_opls()

        for ptype, missing in self.missing.items():
            for miss in missing:
                if miss in self.opls.top[ptype]:
                    self.found.top[ptype][miss] = self.opls.top[ptype][miss].copy()
                elif (miss[::-1] in self.opls.top[ptype] and ptype != 'impropers'):
                    self.found.top[ptype][miss] = self.opls.top[ptype][miss[::-1]].copy()
                else:
                    sys.stderr.write(f"WARNING: Parameter not found ({ptype[0]}) ->  {'  '.join(miss)}\n")
                    self.found.top[ptype][miss] = DynTopol.empty_param(ptype, comment=" WARNING: Parameter not found")

    def write_dynamo(self, file_out:str=None, only_missing:bool=False, ff_file:Union[str,list]=None) -> None:

        if not self.opls: self.gen_opls()
        if only_missing and not self.missing:
            raise ValueError("No missing values provided")
        elif only_missing and not self.found:
            self.find_miss()

        # forcefield = self.found.union(self.opls, 'atomtypes', overwrite=False) if only_missing else self.opls
        forcefield = self.found if only_missing else self.opls

        # new file
        if not ff_file:
            forcefield.union(self.opls, 'residues')
            f = forcefield.build_raw()

        # existing ff topology file
        else:
            ff_external = DynTopol()
            # check single file or list of files
            if isinstance(ff_file, list):
                for i in ff_file:
                    ff_external.union(DynTopol(i), overwrite=True)
            else:
                ff_external.read_ff(ff_file)
            # add residue or modify atom atomtype/charges
            for name, residue in self.opls.top['residues'].items():
                if name not in ff_external.top['residues']:
                    ff_external.top['residues'][name] = residue.copy()
                else:
                    for atom, param in residue['atoms'].items():
                        ff_external.top['residues'][name]['atoms'][atom]['charge'] = param['charge']
                        ff_external.top['residues'][name]['atoms'][atom]['atomtype'] = param['atomtype']
            # add atomtypes / bonds / angles / dihedrals / impropers
            for sect in ('atomtypes', 'bonds', 'angles', 'dihedrals', 'impropers'):
                for name, param in forcefield.top[sect].items():
                    param['comment'] = "Added " + param['comment']
                ff_external.union(forcefield, sect, overwrite=False)
            f = ff_external.build_raw()

        if file_out:
            with open(file_out, 'wt') as outfile:
                outfile.write(f)
        else:
            sys.stdout.write(f)


##  MAIN  #############################################################
if __name__ == '__main__':

    # parser
    parser = argparse.ArgumentParser(prog='itp2dyn',
                                    description='Convert parameters from GROMACS to DYNAMO\n' +
                                                '    .itp (AMBER/GAFF) -> .ff (OPLS)\n',
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('itp', metavar='.itp', type=str, nargs='+',
                        help='input GROMACS parameters files')
    parser.add_argument('-o', metavar='.ff', type=str,
                        help='output DYNAMO parameters file (def: dynamo.ff)')
    parser.add_argument('--ff', metavar='.ff', type=str, nargs='+',
                        help='DYNAMO topology files to take old parameters')
    parser.add_argument("--miss", metavar="", nargs="*",
                        help='reference DYNAMO error message to convert only missing parameters (file or piped)')
    parser.add_argument("--atomtypes", metavar="<>", type=str, default='literal',
                        choices=('literal', 'inter', 'element'),
                        help='mode for atomtypes reading {literal, inter, element} (def: literal)')
    args = parser.parse_args()

    # input file assign
    itp_files = args.itp
    ff_files  = args.ff
    out_file = args.o or 'dynamo.ff'
    missing = [] if args.miss is None else list(fileinput.input(args.miss))
    atomtypes = args.atomtypes

    # parameter conversion
    param = GMXTopology()
    for itp_file in itp_files:
        param.read_itp(itp_file, atomtypes=atomtypes)
    if param.nmolec > 1:
        sys.stderr.write(f"NOTE: {param.nmolec} molecules found -> {' '.join(param.molnames)}\n")
    param.gen_opls()
    if missing: param.read_miss(missing)
    param.write_dynamo(out_file, only_missing=missing, ff_file=ff_files)
