#!/usr/bin/env python3

# Description: Convert parameters from GROMACS' .itp (AMBER/GAFF) to DYNAMO's .ff (OPLS)
# Last update: 26-04-2021

# Inspired from ParmEd (gromacstop.py)
# https://github.com/ParmEd/ParmEd

import argparse
import fileinput
import re
import sys
from copy import deepcopy
from textwrap import dedent
from typing import Union

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
    """
        fDynamo force field topology
    
        Attributes
        ----------
        top : dict of dict

            mm_definitions : str
            electrostatic : float
            lennard_jones : float
            units : str

            atomtypes : dict of dict(str)
                atnum : int
                sigma : float
                epsilon : float
                comment : str
            residues : dict of dict(str)
                natoms : int
                nbonds : int
                nimpropers : int
                atoms : dict of dict(str)
                    atomtype : str
                    charge : float
                bonds : list of tuples
                impropers : list of tuples
            variants : dict of dict(str)
                ndeletes : int
                nadds : int
                ncharges : int
                nbonds : int
                nimpropers : int
                deletes : list
                adds : dict of dict(str)
                    atomtype : str
                    charge : float
                charges : dict(str)
                bonds : list of tuples
                impropers : list of tuples
            links : list of 'variants' dict
            bonds : dict of dict(tuple)
                k : float
                r : float
                comment : str
            angles : dict of dict(tuple)
                k : float
                theta : float
                comment : str
            dihedrals : dict of dict(tuple)
                v : list
                comment : str
            impropers : dict of dict(tuple)
                v : list
                comment : str
    
    """

    statements_def = {'mm_definitions' : "OPLS_AA 1.0",
                      'electrostatics' : 0.5,
                      'lennard_jones'  : 0.5,
                      'units'          : "kcal/mole"}
    statements = tuple(statements_def.keys())
    sections = ('atomtypes', 'residues', 'variants', 'links',
                'bonds', 'angles', 'dihedrals', 'impropers')

    def __init__(self, ff_file:str=None):
        self.ff_file = ""
        # comprehensive topology
        self.top = {key: None for key in self.statements}
        self.top.update({key: dict() for key in self.sections})
        # raw topology (formatted)
        self.raw_top = []
        self.raw_ndx = {key:[0]*2 for key in self.sections}
        # read topology file
        if ff_file is not None:
            self.read_ff(ff_file)

    def __bool__(self) -> bool:
        return self.ff_file or self.raw_top or any([bool(item) for key, item in self.top.items()])

    def __or__(self, other:'DynamoTopology') -> 'DynamoTopology':
        new = deepcopy(self)
        return new.union(other, None)

    def __sub__(self, other:'DynamoTopology') -> 'DynamoTopology':
        new = deepcopy(self)
        return new.difference(other, None)

    def __repr__(self) -> str:
        return self.raw_str

    @property
    def raw_str(self) -> str:
        return "".join(self.raw_top)

    @property
    def residues(self) -> list:
        return list(self.top['residues'].keys())
    
    def union(self, arg:'DynamoTopology', section:str = None) -> 'DynamoTopology':
        """Combine to topology classes"""
        for sect in self.sections:
            if section and section != sect: continue
            self.top[sect].update(arg.top[sect])
        return self

    def difference(self, arg:'DynamoTopology', section:str = None) -> 'DynamoTopology':
        """Substract a topology class"""
        for sect in self.sections:
            if section and section != sect: continue
            for key in arg.top[sect]:
                self.top[sect].pop(key, None)
        return self

    def read_ff(self, ff_file:str, read_format:bool=True) -> None:
        """Read fDynamo topology file (.ff)"""
        self.__init__()
        self.ff_file = ff_file
        # read file as raw
        self._read_ff_raw(ff_file)
        # read file comprehensively
        if read_format:
            self._read_ff_format()

    def write_ff(self, ff_file:str, build_raw:bool=True) -> None:
        """Write fDynamo topology file (.ff)"""
        if build_raw: self.build_raw()
        with open(ff_file, 'wt') as f:
            f.writelines(self.raw_top)

    def _read_ff_raw(self, ff_file:str) -> None:
        """Read fDynamo topology file to raw format (list of lines)"""
        with open(ff_file, 'rt') as f:
            self.raw_top = f.readlines()
            self._update_ndx()

    def _read_ff_format(self, ff_raw:list=None) -> None:
        """Read fDynamo topology comprehensively"""
        ff_raw = ff_raw or self.raw_top  # default internal raw ff if not other specified
        # initializate topology dictionary
        self.top = {key:None for key in self.statements}
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
            comment = "!".join(line[1:]).strip()

            if keyword == 'end': # ------------------------------------
                current_sect = None

            elif keyword == 'mm_definitions': # -----------------------
                self.top['mm_definitions'] = " ".join(content[1:])

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

            elif current_sect == 'residues' or keyword == 'residue':
                # new residue
                if keyword == 'residue':
                    current_sect = 'residues'
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

    def recount_top_n(self) -> None:
        """Recount number of atoms/bonds/impropers for residues and variants/links"""
        for name, top in self.top['residues'].items():
            top['natoms'] = len(top['atoms'])
            top['nbonds'] = len(top['bonds'])
            top['nimpropers'] = len(top['impropers'])
        for name, top in self.top['variants'].items():
            top['ndeletes'] = len(top['deletes'])
            top['nadds'] = len(top['adds'])
            top['ncharges'] = len(top['charges'])
            top['nbonds'] = len(top['bonds'])
            top['nimpropers'] = len(top['impropers'])
        for name, link in self.top['links'].items():
            for top in link:
                top['ndeletes'] = len(top['deletes'])
                top['nadds'] = len(top['adds'])
                top['ncharges'] = len(top['charges'])
                top['nbonds'] = len(top['bonds'])
                top['nimpropers'] = len(top['impropers'])

    def assign_statements_def(self, overwrite:bool=False) -> None:
        """Assign default statements to topology"""
        for i in self.statements:
            self.top[i] = self.statements_def[i] if not self.top[i] or overwrite else self.top[i]

    def build_raw(self) -> str:
        """Generate raw formatted entry from top"""
        self.assign_statements_def(overwrite=False)
        self.recount_top_n()
        raw = ""
        
        # main header
        raw += "!"+"="*79+"\n" + \
               "!"+" "*22+"OPLS MM Definition File for Proteins\n" + \
               "!"+"="*79+"\n\n"
        raw += self.form_field('mm_definitions') + "\n\n"

        
        # atomtypes
        if self.top['atomtypes']:
            raw += "!"+"="*79+"\n" + \
                   "!  Atom Type Definitions\n" + \
                   "!"+"="*79+"\nTypes\n\n"
            for n, name in enumerate(self.top['atomtypes'].keys()):
                raw += self.form_field('atomtypes', name, header_comment=(n==0))
            raw += "\nEnd\n\n"
            # statements
            for sect in ('electrostatics', 'lennard_jones', 'units'):
                raw += self.form_field(sect)
            raw += "\n\n"
        
        # residues / variants / links
        for sect in ('residues', 'variants', 'links'):
            if self.top[sect]:
                raw += "!"+"="*79+"\n" + \
                       "!  {} Definitions\n".format(sect[:-1].capitalize()) + \
                       "!"+"="*79+"\n{}\n\n".format(sect.capitalize())
                for name in self.top[sect].keys():
                    raw += self.form_field(sect, name, header_comment=True)
                raw += "End\n\n"

        # parameters: bonds / angles / dihedrals / impropers
        if self.top['bonds'] or self.top['angles'] or self.top['dihedrals'] or self.top['impropers']:
            raw += "\n!"+"="*79+"\n" + \
                    "!  Parameter Definitions\n" + \
                    "!"+"="*79+"\nParameters\n\n"
            for sect in ('bonds', 'angles', 'dihedrals', 'impropers'):
                if self.top[sect]:
                    raw += "{}\n".format(sect.capitalize())
                    for n, name in enumerate(self.top[sect].keys()):
                        raw += self.form_field(sect, name, header_comment=(n==0))
                    raw += "End\n\n"
            raw += "End\n"

        raw += "End\n"
        self.raw_top = [line+"\n" for line in raw.split("\n")]
        return self.raw_str

    def form_field(self, form:str, key=None, header_comment:bool=False, external=None) -> str:
        """Build a formatted string of a specific field of the topology (statement/section)"""

        if form not in self.statements + self.sections:
            raise ValueError("Unknown field specifier for formatting to")
        elif form in self.sections and not key:
            raise TypeError("Required argument 'key' with 'section' formatting type")

        topology = external or self.top[form]

        # statements
        if form == 'mm_definitions':
            return "MM_Definitions {:<}\n".format(topology)
        elif form == 'electrostatics':
            return "Electrostatics Scale {:<}\n".format(topology)
        elif form == 'lennard_jones':
            return "Lennard_Jones  Scale {:<}\n".format(topology)
        elif form == 'units':
            return "Units {:<}\n".format(topology)

        # sections
        param = topology if external else topology[key]
        if 'comment' in param:
            if param['comment']: param['comment'] = "! " + param['comment']
        if form not in ('residues', 'variants', 'links'): param.update({'key':key})
        if form == 'atomtypes':
            header = "! Atom Name     Atomic Number        Sigma      Epsilon\n" if header_comment else ""
            return header + "{key:<10s}   {atnum:>10}       {sigma:>12.5f} {epsilon:>12.5f}  {comment:<s}\n".format(**param)
        elif form == 'residues':
            text = []
            text.append("Residue {:<s}\n".format(key))
            text.append("{natoms:>4d} {nbonds:>4d} {nimpropers:>4d}\n\n".format(**param))
            for name, atom in param['atoms'].items():
                text.append("{:<4s}  {atomtype:<4s} {charge:>8.4f}\n".format(name, **atom))
            text.append("\n")
            if param['nbonds'] > 0:
                for n, bond in enumerate(param['bonds']):
                    text.append("{:<4s} {:<4s} ; ".format(*bond))
                    if not (n+1) % 6: text.append("\n")
                if (n+1) % 6: text.append("\n")
                text.append("\n")
            if param['nimpropers'] > 0:
                for n, improper in enumerate(param['impropers']):
                    text.append("{:<4s} {:<4s} {:<4s} {:<4s} ; ".format(*improper))
                    if not (n+1) % 3: text.append("\n")
                if (n+1) % 3: text.append("\n")
                text.append("\n")
            if header_comment:
                text.insert(0,"!"+"-"*79+"\n")
                text.insert(2,"!"+"-"*79+"\n")
                text.insert(3, "! # Atoms, bonds and impropers\n")
            return "".join(text)
        elif form == 'variants':
            text = []
            text.append("Variant {:<s}\n".format(key))
            text.append("{ndeletes:>4d} {nadds:>4d} {ncharges:>4d} {nbonds:>4d} {nimpropers:>4d}\n\n".format(**param))
            if param['ndeletes'] > 0:
                for delete in param['deletes']:
                    text.append("{:<4s} ; ".format(delete))
                text.append("\n\n")
            if param['nadds'] > 0:
                for name, add in param['adds'].items():
                    text.append("{:<4s}  {atomtype:<4s} {charge:>8.4f}\n".format(name, **add))
                text.append("\n")
            if param['ncharges'] > 0:
                for name, charge in param['charges'].items():
                    text.append("{:<4s}  {:<}\n".format(name, charge))
                text.append("\n")
            if param['nbonds'] > 0:
                for n, bond in enumerate(param['bonds']):
                    text.append("{:<4s} {:<4s} ; ".format(*bond))
                    if not (n+1) % 6: text.append("\n")
                if (n+1) % 6: text.append("\n")
                text.append("\n")
            if param['nimpropers'] > 0:
                for n, improper in enumerate(param['impropers']):
                    text.append("{:<4s} {:<4s} {:<4s} {:<4s} ; ".format(*improper))
                    if not (n+1) % 3: text.append("\n")
                if (n+1) % 3: text.append("\n")
                text.append("\n")
            if header_comment:
                text.insert(0,"!"+"-"*79+"\n")
                text.insert(2,"!"+"-"*79+"\n")
                text.insert(3, "! # Deletes, adds, charges, bonds and impropers\n")
            return "".join(text)
        elif form == 'links':
            link1 = self.form_field('variants', key, header_comment, param[0]).replace("Variant", "Link")
            link2 = self.form_field('variants', key, False, param[1]).replace(f"Variant {key}", "")
            return link1 + link2
        elif form == 'bonds':
            header = "! Atoms       FC   Equil\n" if header_comment else ""
            return header + "{key[0]:<4s} {key[1]:<4s} {k:>6.2f} {r:>7.3f}  {comment:<s}\n".format(**param)
        elif form == 'angles':
            header = "! Atoms           FC    Equil\n" if header_comment else ""
            return header + "{key[0]:<4s} {key[1]:<4s} {key[2]:<4s} {k:>5.1f} {theta:>8.2f}  {comment:<s}\n".format(**param)
        elif form in ('dihedrals', 'impropers'):
            header = "! Atoms                  V0      V1      V2      V3\n" if header_comment else ""
            return header + "{:<4s} {:<4s} {:<4s} {:<4s} {:>7.3f} {:>7.3f} {:>7.3f} {:>7.3f}  {:>s}\n".format(*param['key'], *param['v'][0:4], param['comment'])

    def raw_prepend(self, section:str, text:str) -> None:
        """Insert a new entry to the raw topology at the beginning of a specific section"""
        self._raw_insert(1, section, text)

    def raw_append(self, section:str, text:str) -> None:
        """Insert a new entry to the raw topology at the end of a specific section"""
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


##  GROMACS TOPOLOGY CLASS  ###########################################

class GMXTopology():

    #TODO: support more parameters types

    def __init__(self, file_inp=None):
        self.atomtypes = dict()
        self.moleculetype = dict()
        self.system = ""
        self.molecules = dict()
        self.opls = DynamoTopology()
        self.missing = dict()
        self.found = DynamoTopology()
        self.notfound = set()
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
            top['residues'][molname] = {'natoms':0, 'nbonds':0, 'nimpropers':0, 'atoms':dict(), 'bonds':[], 'impropers':[]}
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
        self.found = DynamoTopology()
        self.notfound = set()

        if not self.opls: self.gen_opls()

        for ptype, missing in self.missing.items():
            for miss in missing:
                if miss in self.opls.top[ptype]:
                    self.found.top[ptype][miss] = self.opls.top[ptype][miss].copy()
                    self.found.top[ptype][miss]['comment'] += "Added" 
                elif (miss[::-1] in self.opls.top[ptype] and ptype != 'impropers'):
                    self.found.top[ptype][miss] = self.opls.top[ptype][miss[::-1]].copy()
                    self.found.top[ptype][miss]['comment'] += "Added" 
                else:
                    sys.stderr.write(f"WARNING: Parameter not found ({ptype[0]}) ->  {'  '.join(miss)}\n")
                    self.notfound.add(miss)
                    if ptype == 'bonds':
                        self.found.top[ptype][miss] = {'k':0., 'r':0., 'comment':"Added - WARNING: Parameter not found"}
                    elif ptype == 'angles':
                        self.found.top[ptype][miss] = {'k':0., 'theta':0., 'comment':"Added - WARNING: Parameter not found"}
                    elif ptype in ('dihedrals', 'impropers'):
                        self.found.top[ptype][miss] = {'v':[0.]*4, 'theta':0., 'comment':"Added - WARNING: Parameter not found"}

    def write_dynamo(self, file_out:str=None, only_missing:bool=False, ff_file:Union[str,list]=None) -> None:

        if not self.opls: self.gen_opls()
        if only_missing and not self.missing:
            raise ValueError("No missing values provided")
        elif only_missing and not self.found:
            self.find_miss()

        forcefield = self.found if only_missing else self.opls

        # new file
        if not ff_file:
            forcefield.union(self.opls, 'residues')
            f = forcefield.build_raw()

        # existing ff topology file
        else:
            ff_external = DynamoTopology()
            # check single file or list of files
            if isinstance(ff_file, list):
                for i in ff_file:
                    ff_external.union(DynamoTopology(i))
            else:
                ff_external.read_ff(ff_file)
            # add residue or modify atom charges
            for name, residue in self.opls.top['residues'].items():
                if name not in ff_external.top['residues']:
                    ff_external.top['residues'][name] = residue.copy()
                else:
                    for atom, param in residue['atoms'].items():
                        ff_external.top['residues'][name]['atoms'][atom]['charge'] = param['charge']
            # add bonds / angles / dihedrals / impropers
            for sect in ('bonds', 'angles', 'dihedrals', 'impropers'):
                ff_external.union(forcefield, sect)
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
    parser.add_argument('-ff', metavar='.ff', type=str, nargs='+',
                        help='optional DYNAMO topology files to take old parameters')
    parser.add_argument("-miss", metavar="", nargs="*",
                        help='reference DYNAMO error message to convert only missing parameters (file or piped)')
    args = parser.parse_args()

    # input file assign
    itp_files = args.itp
    ff_files  = args.ff
    out_file = args.o or 'dynamo.ff'
    missing = [] if args.miss is None else list(fileinput.input(args.miss))

    # parameter conversion
    param = GMXTopology()
    for itp_file in itp_files:
        param.read_itp(itp_file, elem_simple=missing)
    if param.nmolec > 1:
        sys.stderr.write(f"NOTE: {param.nmolec} molecules found -> {' '.join(param.molnames)}\n")
    param.gen_opls()
    if missing: param.read_miss(missing)
    param.write_dynamo(out_file, only_missing=missing, ff_file=ff_files)
