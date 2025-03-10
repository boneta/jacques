#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
=======================================================================
  DYNAMON configuration parameters for a calculation
=======================================================================

  Class
  -----

    DynnConfig

"""

import math as m
import sys


class DynnConfig:
    """
        DYNAMON configuration class

        Parameters
        ----------
        file : str, optional
            read options from file (.dynn/.jcq)

        Attributes
        ----------
        opt : dict
        atoms : list of list
        constr : list of dict
        selection : dict
            'sele_name' : dict
                'segi' : dict
                    'resi' : list

        Properties
        ----------
        mode : str
            calculation mode
        name : str
            calculation name
        dynn : str
            DYNAMON configuration file name (.dynn)
        out : str
            output file with distances and energies (.out)
        nconstr : int
            number of defined constraints
        npoints : int
            total number of points covered by constraints

        Methods
        -------
        read_file(file)
            read options from DYNAMON formatted file
        read_opt(**kwargs)
            read options from a dictionary / labelled values
        write_dynn(file)
            write formatted options to file (.dynn)
        def_opt(mode, options_set)
            assign default options for a specific mode
        resolve_constr(n_constr, step)
            Resolve the iteration parameters for the first n constraint
        launch()
            launch the calculation to the queue system
    """

    opt_keys = [
        'mode',

        'name', 'out', 'pdb', 'sys', 'bin', 'sele', 'coord', 'ff', 'seq',
        'cores', 'memory',
        'charge', 'multi', 'force_uhf',
        'semiemp', 'gauss', 'func', 'basis',
        'cg_steps', 'cg_tolerance', 'lbfgsb_steps', 'lbfgsb_tolerance',
        'temp', 'md_step', 'equi', 'prod', 'dcd_freq', 'vel',
        'pbc',
        'loc_steps', 'loc_tolerance', 'ts',
        'irc_dir', 'irc_steps', 'irc_dsp',
        'irc_invert', 'irc_dat', 'irc_crd',
        'dcd_stride', 'int_dcd', 'int_wbox', 'int_ions',
        'kie_atomn', 'kie_skip', 'kie_mass', 'kie_hess',

        'exe', 'queue', 'account', 'dry', 'post'
        ]

    constr_keys = [
        'type',

        'symm',
        'atoms',
        'force',
        'dcrd',
        'n',
        'dinit', 'dend', 'step',
        'dist',
        'dfile'
        ]

    # build empty dictionaries
    opt_empty    = dict({key : None for key in opt_keys})
    constr_empty = dict({key : None for key in constr_keys})

    # sets of available options
    opt_keys_up     = {i.upper() for i in opt_keys}
    opt_keys_low    = {i.lower() for i in opt_keys}
    constr_keys_up  = {i.upper() for i in constr_keys}
    constr_keys_low = {i.lower() for i in constr_keys}

    def __init__(self, file=None):
        '''
            Initialization of empty lists and None dictionaries

            Parameters
            ----------
            file : str, optional
                read options from file (.dynn/.jcq)
        '''

        self.opt    = self.opt_empty.copy()
        self.atoms  = []
        self.constr = []
        self.dim    = 0
        self.selection = dict()
        if file is not None:
            self.read_file(file)

    def __str__(self) -> str:
        '''String formatted representation of the configuration file (.dynn)'''
        #FIXME: Make dynn writing more general
        # accumulate to string
        s = ""
        # mode
        if self.mode is not None:
            s += f"{'MODE':<20} {self.mode.upper()}\n\n"
        # properties
        self.opt['name'] = self.name
        self.opt['out']  = self.out
        # general options
        for option in self.opt_keys[1:-5]:
            if self.opt[option] is not None:
                o = self.opt[option]
                o = f"\"{o}\"" if type(o) == str and any(i in o for i in ["/", ","]) else o
                s += f"{option.upper():<20} {o}\n"
        # atoms
        s += "\n"
        for a in self.atoms:
            s += f"{'ATOM':<20} {a}\n"
        # constraints
        for c in self.constr:
            s += f"\n{'CONSTR':<20} {self._swap_constrtype(c['type'])}\n"
            for option in self.constr_keys[1:]:
                if c[option] is not None:
                    o = c[option]
                    o = f"\"{o}\"" if type(o) == str and "/" in o else o
                    s += f"  {option.upper():<18} {o}\n"
            s += "C\n"
        # selections
        for sele_name, sele in self.selection.items():
            s += f"\nSELECTION {sele_name}\n"
            for segi, resis in sele.items():
                s += " "*4+f"S {segi}\n"
                for resi, atoms in resis.items():
                    s += " "*8+f"R {resi}\n"
                    for name in atoms:
                        s += " "*12+f"A {name}\n"
            s += "SELECTION\n"
        # return final string
        return s

    @property
    def mode(self):
        '''Calculation mode'''
        if self.opt['mode'] is not None:
            return self.opt['mode'].lower()
        else:
            return None

    @property
    def name(self):
        '''Calculation name'''
        if self.opt['name'] is not None:
            return self.opt['name']
        elif self.opt['mode'] is not None:
            return self.mode.upper()
        else:
            return None

    @property
    def dynn(self):
        '''DYNAMON configuration file name (.dynn)'''
        return self.name + '.dynn'

    @property
    def out(self):
        '''Output file with distances and energies (.out)'''
        if self.opt['out'] is not None:
            return self.opt['out']
        else:
            return self.name + ".out"

    @property
    def nconstr(self):
        '''Number of defined constraints'''
        return len(self.constr)

    @property
    def npoints(self):
        '''Total number of points covered by constraints'''
        if self.dim == 0:
            return 0
        else:
            n = 1
            for d in range(self.dim):
                n *= self.constr[d]['n']+1
            return n

    def read_file(self, file):                  #FIXME: Make reading more robust
        '''
            Read options file in DYNAMON format

            Parameters
            ----------
            file : str
                options file to be read
        '''

        # read all file
        with open(file, 'rt') as inputfile:
            opts = inputfile.readlines()
            # remove comment or empty lines and space-split the first word
            opts = [list(map(str.strip, line.strip().split(' ',1))) for line in opts
                    if line.strip() and not line.startswith(("!","#"))]

        # loop line by line to assign options
        n = 0
        while n < len(opts):
            line = list(map(str.strip, opts[n]))
            if line[0] in (self.opt_keys_up|self.opt_keys_low):
                self.opt[line[0].lower()] = line[1].strip('"')
            elif line[0].upper() in ('ATOM', 'A'):
                self.atoms.append(line[1])
            elif line[0].upper() in ('CONSTR', 'C'):
                self.constr.append(self.constr_empty.copy())
                self.constr[-1]['type'] = self._swap_constrtype(line[1])      # type of constr
                # read constraint options
                while n+1 < len(opts):
                    n += 1
                    line = opts[n]
                    if line[0] == 'CONSTR' or line[0] == 'C': break
                    if line[0] in (self.constr_keys_up|self.constr_keys_low):
                        self.constr[-1][line[0].lower()] = line[1].strip('"')
                    else:
                        print("WARNING: Unrecognized option:", line[0])
            elif line[0].upper() == 'SELECTION':
                sele_name = line[1].split()[0].upper()
                # read whole selection to dict of dict of list
                sele = dict()
                segi = ""
                resi = ""
                while n+1 < len(opts):
                    n += 1
                    line = opts[n]
                    if line[0].upper() == 'SELECTION': break
                    subsect = line[0].upper()
                    select = line[1].upper()
                    if subsect == "S":
                        segi = select
                        resi = ""
                        sele[select] = dict()
                    elif subsect == "R":
                        resi = int(select)
                        sele[segi][int(select)] = []
                    elif subsect == "A":
                        sele[segi][resi].append(select)
                self.selection[sele_name] = sele
            else:
                print("WARNING: Unrecognized option:", line[0])
            n += 1

    def read_opt(self, **kwargs):
        '''Read options from a dictionary / labelled values'''
        # read configuration file first
        if 'f' in kwargs and kwargs['f'] is not None:
            self.read_file(kwargs['f'])
        for option, value in kwargs.items():
            if option == 'dim':
                self.dim = value
            elif option in self.opt_keys:
                self.opt[option] = value or self.opt[option]

    def write_dynn(self, file=None):
        '''
            Write formatted options to file (.dynn)

            Parameters
            ----------
            file : str, optional
                file name to be written
                if None, the default is used based on name/mode
        '''
        file = file or self.dynn
        with open(file, 'wt') as f:
            f.write(str(self))

    def def_opt(self, mode, opt_settings):
        '''
            Assign default options for a specific mode

            Parameters
            ----------
            mode : str
                calculation mode
            opt_settings : dic
                dictionary of options in the format of 'settings.yaml'
        '''

        if mode is None:
            sys.stdout.write("WARNING: Default settings not assigned. Missing MODE\n")
            return

        # general options of specific mode and general ("all")
        for def_dic in [mode, "all"]:
            if def_dic not in opt_settings.keys():
                sys.stdout.write("WARNING: Default settings not found for '{}'\n".format(def_dic))
                continue
            else:
                for option in opt_settings[def_dic].keys():
                    if option in self.opt_keys[1:] and self.opt[option] is None:
                        self.opt[option] = opt_settings[def_dic][option]

            # constraints
            if 'constr' in opt_settings[def_dic].keys():
                for c in self.constr:
                    for option in opt_settings[def_dic]['constr'].keys():
                        if option in self.constr_keys[1:] and c[option] is None:
                            c[option] = opt_settings[def_dic]['constr'][option]
        # dimensions
        if self.dim == 0:
            self.dim = self.nconstr

    def resolve_constr(self, n_constr=None, step=0.05):
        '''
            Resolve the iteration parameters for the
            first nth constraint: n, dinit, step & dend

            At least three are needed to resolve the fourth,
            altough 'step' can be taken from default.
            In case of conflict, 'step' is overwritten.

            If 'dcrd' is set, that constraint is skipped.

            Parameters
            ----------
            n_constr : int, optional
                number of constraints to resolve, default 'None' for all
            step : float, optional
                default step lenght if not found (def: 0.05)
        '''

        # check number of constraints
        n_constr = n_constr or self.nconstr
        if n_constr > self.nconstr:
            raise ValueError(f"More contraints requested to resolve ({n_constr}) than defined ({self.nconstr})")

        types_dict = {'dinit':float, 'dend':float, 'step':float, 'n':int}

        for nth, c in enumerate(self.constr[0:n_constr]):
            # type conversion
            for i,j in types_dict.items():
                if c[i] is not None: c[i] = j(c[i])
            # default step
            if not c['step']: c['step'] = step
            # check enough defined parameters to solve
            if c['dcrd']:
                continue
            elif (c['dinit'], c['dend'], c['n']).count(None) > 1:
                raise ValueError(f"Not enough parameters specified to resolve constraint {nth+1}")
            # dinit & dend -> n [& step]
            if c['dinit'] and c['dend']:
                diff = c['dend'] - c['dinit']
                if not c['n']:
                    c['n'] = m.ceil(diff/c['step'])
                    if c['n'] < 0: c['n'], c['step'] = -c['n'], -c['step']
                else:
                    c['step'] = round(diff/c['n'], 3)
            # dinit & n [& step] -> dend
            elif c['dinit']:
                c['dend'] = c['dinit'] + c['step']*c['n']
            # dend & n [& step] -> dinit
            elif c['dend']:
                c['dinit'] = c['dend'] - c['step']*c['n']

    @staticmethod
    def _swap_constrtype(c):
        '''
            Convert between constraint notation types: str <-> int

            Parameters
            ----------
            c : str or int
                constraint type to convert

            Returns
            -------
            y : int or str
                converted constraint of the opposite type
        '''

        str2int = { 'D':1, 'M':2 }
        int2str = { 1:'D', 2:'M' }

        # string to integer
        if isinstance(c,str):
            return str2int[c.upper()]
        # integer to string
        elif isinstance(c,int):
            return int2str[c]

    def launch(self):
        '''Launch the calculation'''
        from jacques.calculation_modes import launch_calc
        launch_calc(self)
