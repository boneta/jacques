#!/usr/bin/env python3

# Description: Convert parameters from GROMACS' .itp (AMBER/GAFF) to DYNAMO's .ff (OPLS)
# Last update: 22-04-2021

# Inspired from ParmEd (gromacstop.py)
# https://github.com/ParmEd/ParmEd

import argparse
import fileinput
import re
import sys
from textwrap import dedent

from amber2opls import DihedralParam


##  INFO  #############################################################

"""
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
ptable_simple_inv = { n:elem for elem, n in ptable_simple.items() }


##  DYNAMO TOPOLOGY CLASS  ############################################

class DynamoTopology():
    """fDynamo force field topology"""

    sections = ['atomtypes', 'residues', 'variants', 'links',
                'bonds', 'angles', 'dihedrals', 'impropers']

    def __init__(self, ff_file:str=None):
        self.ff_file = ""
        self.top = {key:dict() for key in self.sections}
        self.raw_top = []
        self.raw_ndx = {key:[0]*2 for key in self.sections}
        # read topology file
        if ff_file is not None:
            self.read_ff(ff_file)

    def __bool__(self):
        return bool(self.ff_file)

    def __repr__(self):
        return self.raw_str

    @property
    def raw_str(self):
        return "".join(self.raw_top)

    def read_ff(self, ff_file:str, read_format:bool=True) -> None:
        """Read fDynamo topology file (.ff)"""
        self.__init__()
        self.ff_file = ff_file
        # read file as raw
        self._read_ff_raw(ff_file)
        # read file comprehensively
        if read_format:
            self._read_ff_format()

    def _read_ff_raw(self, ff_file:str) -> None:
        """Read fDynamo topology file to raw format (list of lines)"""
        with open(ff_file, 'rt') as f:
            self.raw_top = f.readlines()
            self._update_ndx()

    def _read_ff_format(self, ff_raw:list=None) -> None:
        """Read fDynamo topology comprehensively"""
        # default internal raw ff if not other specified
        ff_raw = ff_raw if ff_raw is not None else self.raw_top
        # initializate topology dictionary
        self.top = {'ff':None, 'mm_definitions':None, 'electrostatics':None, 'lennard_jones':None, 'units':None}
        self.top.update({key:dict() for key in self.sections})
        # remove comment lines and blank lines
        topology = [line.strip() for line in ff_raw if line.strip() and not line.startswith("!")]
        # read line by line and assign depending on the section
        # based on fDynamo method: MM_FILE_PROCESS() @ mm_file_io.F90
        current_sect = None
        for line in topology:
            line = line.split("!")
            keyword = line[0].split()[0].lower()
            content = line[0].split()
            comment = "!".join(line[1:])

            if keyword == 'end': # ------------------------------------
                current_sect = None

            elif keyword == 'mm_definitions': # -----------------------
                self.top['ff'] = content[1:]

            elif keyword in ('electrostatics', 'lennard_jones'): # ----
                if content[1].lower() != "scale":
                    sys.stderr.write(f"WARNING: Unrecognized {keyword.upper()} keyword\n")
                self.top[keyword] = float(content[2])

            elif keyword == 'units': # --------------------------------
                self.top['units'] = str(content[1])

            elif keyword in self.sections or keyword == 'types': # ----
                current_sect = keyword

            elif current_sect == 'types': # ---------------------------
                attype = str(content[0]).upper()
                param = { 'atnum' : int(content[1]),
                          'sigma' : float(content[2]),
                          'epsilon' : float(content[3]),
                          'comment' : comment }
                self.top['atomtypes'][attype] = param

            elif current_sect == 'residues': # ------------------------
                # new residue
                if keyword == 'residue':
                    name = content[1].upper()
                    self.top['residues'][name] = dict()
                    continue
                # read residue properties
                top = self.top['residues'][name]
                if not top:
                    top['natoms'] = int(content[0])
                    top['nbonds'] = int(content[1])
                    top['nimpropers'] = int(content[2])
                    top['atoms'] = dict()
                    top['bonds'] = []
                    top['impropers'] = []
                elif len(top['atoms']) < top['natoms']:
                    param = { 'atomtype' : str(content[1]).upper(),
                              'charge' : float(content[2]) }
                    top['atoms'][str(content[0]).upper()] = param
                elif len(top['bonds']) < top['nbonds']:
                    bonds = [(b.split()[0].upper(), b.split()[1].upper())
                              for b in " ".join(content).split(";") if b.strip()]
                    top['bonds'].extend(bonds)
                elif len(top['impropers']) < top['nimpropers']:
                    impropers = [(i.split()[0].upper(), i.split()[1].upper(),
                                  i.split()[2].upper(), i.split()[3].upper())
                                  for i in " ".join(content).split(";") if i.strip()]
                    top['impropers'].extend(impropers)

            elif current_sect in ('variants', 'links'): # -------------
                # new variant
                if keyword == 'variant':
                    name = " ".join(content[1:]).upper()
                    self.top['variants'][name] = dict()
                    continue
                # new link
                elif keyword == 'link':
                    name = " ".join(content[1:]).upper()
                    self.top['links'][name] = [dict()]*2
                    continue
                # read variant/link properties
                if current_sect == 'variants':
                    top = self.top['variants'][name]
                elif current_sect == 'links':
                    top = self.top['links'][name][0]
                    fields = ('deletes','adds','charges','bonds','impropers')
                    # check if first link is complete
                    if top and all([ top['n'+field]==len(top[field]) for field in fields ]):
                        top = self.top['links'][name][1]
                if not top:
                    top['ndeletes'] = int(content[0])
                    top['nadds'] = int(content[1])
                    top['ncharges'] = int(content[2])
                    top['nbonds'] = int(content[3])
                    top['nimpropers'] = int(content[4])
                    top['deletes'] = []
                    top['adds'] = dict()
                    top['charges'] = dict()
                    top['bonds'] = []
                    top['impropers'] = []
                elif len(top['deletes']) < top['ndeletes']:
                    deletes = [b.strip().upper()
                               for b in " ".join(content).split(";") if b.strip()]
                    top['deletes'].extend(deletes)
                elif len(top['adds']) < top['nadds']:
                    param = { 'atomtype' : str(content[1]).upper(),
                              'charge' : float(content[2]) }
                    top['adds'][str(content[0]).upper()] = param
                elif len(top['charges']) < top['ncharges']:
                    top['charges'][str(content[0]).upper()] = float(content[1])
                elif len(top['bonds']) < top['nbonds']:
                    bonds = [(b.split()[0].upper(), b.split()[1].upper())
                              for b in " ".join(content).split(";") if b.strip()]
                    top['bonds'].extend(bonds)
                elif len(top['impropers']) < top['nimpropers']:
                    impropers = [(i.split()[0].upper(), i.split()[1].upper(),
                                  i.split()[2].upper(), i.split()[3].upper())
                                  for i in " ".join(content).split(";") if i.strip()]
                    top['impropers'].extend(impropers)

            elif keyword == 'parameters': # ---------------------------
                continue

            elif current_sect == 'bonds': # ---------------------------
                atoms = tuple(str(a).upper() for a in content[:2])
                param = { 'k' : float(content[2]),
                          'r' : float(content[3]),
                          'comment' : comment }
                self.top['bonds'][atoms] = param

            elif current_sect == 'angles': # --------------------------
                atoms = tuple(str(a).upper() for a in content[:3])
                param = { 'k' : float(content[3]),
                          'theta' : float(content[4]),
                          'comment' : comment }
                self.top['angles'][atoms] = param

            elif current_sect in ('dihedrals', 'impropers'): # --------
                atoms = tuple(str(a).upper() for a in content[:4])
                param = { 'v' : [float(v) for v in content[4:]],
                          'comment' : comment }
                self.top[current_sect][atoms] = param

            else:
                sys.stderr.write("WARNING: Unrecognized option in FF file -> {}\n".format("!".join(line)))

    def _update_ndx(self) -> None:
        """Update line indexes of raw stored topology"""
        current_sect = None
        for n, line in enumerate(self.raw_top):
            # remove comments and empty lines
            if not line.strip() or line.startswith("!"):
                continue
            keyword = line.split()[0].lower()
            if keyword in self.sections or keyword == 'types':
                current_sect = keyword if keyword != 'types' else 'atomtypes'
                self.raw_ndx[current_sect][0] = n
            elif keyword == 'end' and current_sect is not None:
                self.raw_ndx[current_sect][1] = n
                current_sect = None

    def raw_prepend(self, section:str, text:str) -> None:
        self._raw_insert(1, section, text)

    def raw_append(self, section:str, text:str) -> None:
        self._raw_insert(-1, section, text)

    def _raw_insert(self, action:int, section:str, text:str) -> None:
        """Insert a new entry to the raw topology in a specific section"""
        if section not in self.sections or abs(action) != 1:
            raise ValueError
        ndx = self.raw_ndx[section]
        for line in text.splitlines()[::-1]:
            line += "\n"
            # prepend
            if action == 1:
                    self.raw_top.insert(ndx[0]+1, line)
            # append
            elif action == -1:
                self.raw_top.insert(ndx[1], line)
        self._update_ndx()


##  FORCE FIELS PARAMTERS CLASS  ######################################

class ffParameters():

    #TODO: support more parameters types

    def __init__(self, file_inp=None):
        self.atomtypes = dict()
        self.moleculetype = dict()
        self.system = ""
        self.molecules = dict()
        self.opls = dict()
        self.missing = dict()
        self.found = dict()
        self.notfound = []
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

    def read_itp(self, file_inp:str, elem_simple:bool=False) -> None:

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

        # process line by line
        attype_re = re.compile(r'[a-zA-Z]') if elem_simple else re.compile(r'.+')
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
                attype = attype_re.findall(str(words[-1]))[0]
                param['epsilon'] = float(words[0])
                param['sigma'] = float(words[1])
                param['ptype'] = str(words[2])
                param['charge'] = float(words[3])
                param['mass'] = float(words[4])
                # guess atomic number if absent
                if not param['atnum']: param['atnum'] = ptable_simple[attype[0].upper()]
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
                param['attype'] = attype_re.findall(str(words[1]))[0]
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

        # dictionaries -> if a parameters repeats its is overwritten
        self.opls = { 'atomtypes' : dict(),
                      'residues' : dict(),
                      'bonds' : dict(),
                      'angles' : dict(),
                      'dihedrals' : dict(),
                      'impropers' : dict() }

        # atom types
        for attype, atom in self.atomtypes.items():
            self.opls['atomtypes'][attype.upper()] = [atom['atnum'], atom['sigma']*10, atom['epsilon']]

        # atoms per residue and unified bonds/angles/dihedrals/impropers
        param = DihedralParam()
        for molname, molec in self.moleculetype.items():
            # atoms in residue
            self.opls['residues'][molname] = dict()
            for nr, atom in molec['atoms'].items():
                if atom['name'] in self.opls['residues'][molname]:
                    raise ValueError(f"opls conversion -> duplicated atom definition (residue {molname}, atom {atom['name']})")
                self.opls['residues'][molname][atom['name']] = [atom['attype'].upper(), atom['charge']]
            # bonds
            for atoms, bond in molec['bonds'].items():
                attypes = tuple(molec['atoms'][nr]['attype'].upper() for nr in atoms)
                self.opls['bonds'][attypes] = [bond['k']/1000, bond['r']*10]
            # angles
            for atoms, angle in molec['angles'].items():
                attypes = tuple(molec['atoms'][nr]['attype'].upper() for nr in atoms)
                self.opls['angles'][attypes] = [angle['cth']/10, angle['theta']]
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
                    self.opls[param_type][attypes] = param.p['opls']

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
        self.found = { 'bonds' : dict(),
                       'angles' : dict(),
                       'dihedrals' : dict(),
                       'impropers' : dict() }
        self.notfound = set()

        if not self.opls: self.gen_opls()

        for ptype, missing in self.missing.items():
            for miss in missing:
                if miss in self.opls[ptype]:
                    self.found[ptype][miss] = self.opls[ptype][miss].copy()
                elif (miss[::-1] in self.opls[ptype] and ptype != 'impropers'):
                    self.found[ptype][miss] = self.opls[ptype][miss[::-1]].copy()
                else:
                    sys.stderr.write(f"WARNING: Parameter not found ({ptype[0]}) ->  {'  '.join(miss)}\n")
                    n_zeros = 2 if ptype in ('bonds', 'angles') else 5 if ptype in ('dihedrals', 'impropers') else 0
                    self.found[ptype][miss] = [0]*n_zeros
                    self.notfound.add(miss)

    def write_dynamo(self, file_out:str=None, only_missing:bool=False, ff_file:str=None) -> None:

        #TODO: less redundant method between new file and append

        def _check_notfound(atoms):
            return " "*4+"! WARNING: Parameter not found" if atoms in self.notfound else ""

        if not self.opls:
            self.gen_opls()
        if only_missing and not self.missing:
            raise ValueError("No missing values provided")
        elif only_missing and not self.found:
            self.find_miss()

        # NEW FILE ----------------------------------------------------
        if not ff_file:
            f = []
            f.extend("!"+"="*79+"\n"+
                    "!"+" "*22+"OPLS MM Definition File for Proteins\n"+
                    "!"+"="*79+"\n"+
                    "!"+" "*22+f"Automatically generated from: {self.file_inp}\n"+
                    "!"+"-"*79+"\n\n")
            f.extend(dedent("""
                                MM_Definitions OPLS_AA 1.0\n
                                Electrostatics Scale 0.5
                                Lennard_Jones  Scale 0.5\n
                                Units kcal/mole\n\n
                            """))

            # atomtypes
            if self.opls['atomtypes'] and not only_missing:
                f.extend("!"+"="*79+"\n" +
                        "!  Atom Type Definitions\n" +
                        "!"+"="*79+"\n"+
                        "Types\n\n"+
                        "! Atom Name     Atomic Number        Sigma      Epsilon\n")
                for attype, atom in self.opls['atomtypes'].items():
                    f.extend("{:<10s}   {:>10}       {:>12.5f} {:>12.5f}\n".format(attype, *atom))
                f.extend("\nEnd\n\n")

            # residues
            if self.opls['residues']:
                f.extend("!"+"="*79+"\n" +
                        "!  Residue Definitions\n" +
                        "!"+"="*79+"\n" +
                        "Residues\n")
                for resname, atoms in self.opls['residues'].items():
                    f.extend("\n!"+"-"*79+"\n" +
                            "Residue {}\n".format(resname) +
                            "!"+"-"*79+"\n" +
                            "! # Atoms, bonds and impropers\n" +
                            "{:>3d}{:>3d}{:>3d}\n".format(len(atoms),0,0))
                    for name, atom in atoms.items():
                        f.extend("{:<4s}  {:<4s} {:>6.2f}\n".format(name, *atom))
                f.extend("\nEnd\n")

            if only_missing:
                self.opls = self.found

            # parameters
            f.extend("\n\n"+
                    "!"+"="*79+"\n" +
                    "!  Parameter Definitions\n" +
                    "!"+"="*79+"\n" +
                    "Parameters\n")

            # bonds
            if self.opls['bonds']:
                f.extend("\nBonds\n"+
                        "! Atoms       FC   Equil\n")
                for atoms, bond in self.opls['bonds'].items():
                    f.extend("{:<4s} {:<4s} {:>6.2f} {:>7.3f} {:<10s}\n".format(
                            *atoms, *bond, _check_notfound(atoms)))
                f.extend("End\n")
            # angles
            if self.opls['angles']:
                f.extend("\nAngles\n"+
                        "! Atoms           FC    Equil\n")
                for atoms, angle in self.opls['angles'].items():
                    f.extend("{:<4s} {:<4s} {:<4s} {:>5.1f}  {:>8.2f} {:<10s}\n".format(
                            *atoms, *angle, _check_notfound(atoms)))
                f.extend("End\n")
            # dihedrals
            if self.opls['dihedrals']:
                f.extend("\nDihedrals\n"+
                        "! Atoms                  V0      V1      V2      V3\n")
                for atoms, dihedral in self.opls['dihedrals'].items():
                    f.extend("{:<4s} {:<4s} {:<4s} {:<4s} {:>7.3f} {:>7.3f} {:>7.3f} {:>7.3f} {:>10s}\n".format(
                            *atoms, *dihedral[:-1], _check_notfound(atoms)))
                f.extend("End\n")
            # impropers
            if self.opls['impropers']:
                f.extend("\nImpropers\n"+
                        "! Atoms                  V0      V1      V2      V3\n")
                for atoms, improper in self.opls['impropers'].items():
                    f.extend("{:<4s} {:<4s} {:<4s} {:<4s} {:>7.3f} {:>7.3f} {:>7.3f} {:>7.3f} {:>10s}\n".format(
                            *atoms, *improper[:-1], _check_notfound(atoms)))
                f.extend("End\n")

            f.extend("\nEnd\nEnd\n")
            f = "".join(f)

        # EXISTING FF TOPOLOGY ----------------------------------------
        else:
            ff = DynamoTopology(ff_file)
            # residues
            for resname, atoms in reversed(self.opls['residues'].items()):
                f = []
                f.extend("\n!"+"-"*79+"\n" +
                         "Residue {}\n".format(resname) +
                         "!"+"-"*79+"\n" +
                         "! # Atoms, bonds and impropers\n" +
                         "{:>3d}{:>3d}{:>3d}\n".format(len(atoms), 0, 0))
                for name, atom in atoms.items():
                    f.extend("{:<4s}  {:<4s} {:>6.2f}\n".format(name, *atom))
                ff.raw_prepend('residues', "".join(f))

            if only_missing:
                self.opls = self.found

            # bonds
            f = ["! Additional parameters\n"]
            for atoms, bond in self.opls['bonds'].items():
                f.extend("{:<4s} {:<4s} {:>6.2f} {:>7.3f} {:<10s}\n".format(
                         *atoms, *bond, _check_notfound(atoms)))
            ff.raw_append('bonds', "".join(f))

            # angles
            f = ["! Additional parameters\n"]
            for atoms, angle in self.opls['angles'].items():
                f.extend("{:<4s} {:<4s} {:<4s} {:>5.1f}  {:>8.2f} {:<10s}\n".format(
                         *atoms, *angle, _check_notfound(atoms)))
            ff.raw_append('angles', "".join(f))

            # dihedrals
            f = ["! Additional parameters\n"]
            for atoms, dihedral in self.opls['dihedrals'].items():
                f.extend("{:<4s} {:<4s} {:<4s} {:<4s} {:>7.3f} {:>7.3f} {:>7.3f} {:>7.3f} {:>10s}\n".format(
                         *atoms, *dihedral[:-1], _check_notfound(atoms)))
            ff.raw_append('dihedrals', "".join(f))

            # impropers
            f = ["! Additional parameters\n"]
            for atoms, improper in self.opls['impropers'].items():
                f.extend("{:<4s} {:<4s} {:<4s} {:<4s} {:>7.3f} {:>7.3f} {:>7.3f} {:>7.3f} {:>10s}\n".format(
                         *atoms, *improper[:-1], _check_notfound(atoms)))
            ff.raw_append('impropers', "".join(f))

            f = ff.raw_str

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
                        help='input GROMACS parameters file')
    parser.add_argument('-o', metavar='.ff', type=str,
                        help='output DYNAMO parameters file (def: dynamo.ff)')
    parser.add_argument('-ff', metavar='.ff', type=str,
                        help='optional DYNAMO topology file to include the new parameters')
    parser.add_argument("-miss", metavar="", nargs="*",
                        help='reference DYNAMO error message to convert only missing parameters (file or piped)')
    args = parser.parse_args()

    # input file assign
    itp_files = args.itp
    ff_file  = args.ff
    out_file = args.o if args.o is not None else 'dynamo.ff'
    missing = [] if args.miss is None else list(fileinput.input(args.miss))

    # parameter conversion
    param = ffParameters()
    for itp_file in itp_files:
        param.read_itp(itp_file, elem_simple=missing)
    if param.nmolec > 1:
        sys.stderr.write(f"NOTE: {param.nmolec} molecules found -> {' '.join(param.molnames)}\n")
    param.gen_opls()
    if missing: param.read_miss(missing)
    param.write_dynamo(out_file, only_missing=missing, ff_file=ff_file)
