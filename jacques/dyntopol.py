#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
=======================================================================
  fDynamo topology handling
=======================================================================

  Class
  -----

    DynTopol

"""

import sys
from copy import deepcopy


class DynTopol():
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
        # missing parameters
        self.missing = {'atomtypes':{}, 'bonds':{}, 'angles':{}, 'dihedrals':{}, 'impropers':{}}
        # raw topology (formatted)
        self.raw_top = []
        self.raw_ndx = {key:[0]*2 for key in self.sections}
        # read topology file
        if ff_file is not None:
            self.read_ff(ff_file)

    def __bool__(self) -> bool:
        return self.ff_file or self.raw_top or any([bool(item) for key, item in self.top.items()])

    def __or__(self, other:'DynTopol') -> 'DynTopol':
        new = deepcopy(self)
        return new.union(other, None, True)

    def __sub__(self, other:'DynTopol') -> 'DynTopol':
        new = deepcopy(self)
        return new.difference(other, None, True)

    def __repr__(self) -> str:
        return self.raw_str

    @property
    def raw_str(self) -> str:
        return "".join(self.raw_top)

    @property
    def residues(self) -> list:
        return list(self.top['residues'].keys())

    def union(self, other:'DynTopol', section:str=None, overwrite:bool=True) -> 'DynTopol':
        """Combine to topology classes"""
        for sect in self.sections:
            if section and section != sect: continue
            if overwrite:
                self.top[sect].update(other.top[sect])
            else:
                for name, param in other.top[sect].items():
                    self.top[sect][name] = self.top[sect].get(name, param)
        self.remove_duplicates(section, overwrite)
        return self

    def difference(self, other:'DynTopol', section:str=None) -> 'DynTopol':
        """Substract a topology class"""
        for sect in self.sections:
            if section and section != sect: continue
            for key in other.top[sect]:
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

    def remove_duplicates(self, section:str=None, overwrite:bool=True) -> None:
        """Remove duplicate entries of bonds/angles/dihedrals due to inverted definitions"""
        for sect in ('bonds', 'angles', 'dihedrals'):
            if section and section != sect: continue
            keys = list(self.top[sect].keys())
            n = 0
            while n < len(keys):
                key = keys[n]
                yek = tuple(reversed(key))
                n += 1
                if yek == key:  # skip palindromic
                    continue
                elif yek in keys:
                    if overwrite: self.top[sect][key] = self.top[sect][yek]
                    self.top[sect].pop(yek)
                    keys.remove(yek)

    def find_missing(self) -> dict:
        """Find missing paramters from residues (atomtypes/bonds/impropers)"""
        self.missing = {'atomtypes':{}, 'bonds':{}, 'angles':{}, 'dihedrals':{}, 'impropers':{}}
        # atomtypes
        atomtypes_res = {p_atom['atomtype'] for resname, p_res in self.top['residues'].items() for atom, p_atom in p_res['atoms'].items()}
        self.missing['atomtypes'] = atomtypes_res - set(self.top['atomtypes'].keys())
        # bonds / dihedrals
        for sect in ('bonds', 'impropers'):
            param_res = {name for resname, p_res in self.top['residues'].items() for name in p_res[sect]}
            param_top = set(self.top[sect].keys())
            if sect != 'impropers': param_top.update([name[::-1] for name in self.top[sect].keys()])
            self.missing[sect] = param_res - param_top
        return self.missing

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

    @staticmethod
    def empty_param(section:str, **kwargs) -> dict:
        """Return a empty parameter of a specified section"""
        if section == 'atomtypes':
            param = {'atnum':0, 'sigma':0.0, 'epsilon':0.0, 'comment':""}
        elif section == 'residues':
            param = {'natoms':0, 'nbonds':0, 'nimpropers':0, 'atoms':dict(), 'bonds':[], 'impropers':[]}
        elif section in ('variants', 'links'):
            param = {'ndeletes':0, 'nadds':0, 'ncharges':0, 'nbonds':0, 'nimpropers':0,
                     'deletes':[], 'adds':dict(), 'charges':[], 'bonds':[], 'impropers':[]}
        elif section == 'bonds':
            param = {'k':0.0, 'r':0.0, 'comment':""}
        elif section == 'angles':
            param = {'k':0.0, 'theta':0.0, 'comment':""}
        elif section in ('dihedrals', 'impropers'):
            param = {'v':[0.0]*5, 'comment':""}
        else:
            raise ValueError("Unknown section")
        param.update(kwargs)    # assign known arguments
        return param
