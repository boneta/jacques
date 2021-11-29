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

import os
import sys
import re
import glob
import shutil
import itertools
import math as m
from textwrap import dedent

try:
    import jacques.queues as queues
    jacques_import = True
except ImportError:
    jacques_import = False

class DynnConfig:
    '''
        DYNAMON configuration class

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
        nconstr : int
            number of defined constraints
        npoints : int
            total number of points covered by constraints

        Methods
        -------
        __init__(file)
            Initialization of empty lists and None dictionaries
        read_file(file)
            read options from DYNAMON formatted file
        write(file, constr)
            write formatted options to screen or file
        def_opt(mode, options_set)
            assign default options for a specific mode
        resolve_constr(n_constr, step)
            Resolve the iteration parameters for the first n constraint
        guess_exe_bin(dynamon_path)
            guess either executable or binary from the other
        launch()
            launch the calculation to the queue system
        post(rm)
            post-process routine after a calculation
    '''

    opt_keys = [
        'mode',

        'name', 'out', 'sys', 'bin', 'sele', 'coord', 'ff', 'seq',
        'cores', 'memory',
        'charge', 'multi', 'force_uhf',
        'semiemp', 'gauss', 'func', 'basis',
        'cg_steps', 'cg_tolerance', 'lbfgsb_steps', 'lbfgsb_tolerance',
        'temp', 'md_step', 'equi', 'prod', 'dcd_freq', 'vel',
        'pbc',
        'loc_steps', 'loc_tolerance', 'ts',
        'irc_dir', 'irc_steps', 'irc_dsp',
        'dcd_stride', 'int_dcd', 'int_wbox', 'int_ions',
        'kie_atomn', 'kie_skip', 'kie_mass', 'kie_hess',

        'exe', 'queue', 'dry', 'post'
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
                read options from file
        '''

        self.opt    = self.opt_empty.copy()
        self.atoms  = []
        self.constr = []
        self.dim    = 0
        self.selection = dict()
        if file is not None:
            self.read_file(file)

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
        else:
            return self.mode.upper()

    @property
    def out(self):
        '''Output file with distances and energies'''
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
        '''
            Read options from a dictionary / labelled values
        '''

        for option, value in kwargs.items():
            if not value:
                continue
            elif option in ('f', 'file'):
                self.read_file(value)
            elif option == 'dim':
                self.dim = value
            elif option in self.opt_keys:
                self.opt[option] = value

    def write(self, file=None, constr=True):   #FIXME: Make writing more general
        '''
            Write formatted options

            Parameters
            ----------
            file : str
                file write (.dynn)
                if None, print to the screen
            constr : logical
                write constraints secctions
        '''

        # accumulate to string and write at the end
        s = ""

        # mode
        if self.mode is not None:
            s += "{:<20} {}\n\n".format("MODE", self.mode.upper())
        else:
            sys.stdout.write("WARNING: 'MODE' parameter not specified\n")

        # general options
        for option in self.opt_keys[1:-4]:
            if self.opt[option] is not None:
                if type(self.opt[option])==str and ("/" in self.opt[option] or "," in self.opt[option]):
                    s += "{:<20} \"{}\"\n".format(option.upper(),self.opt[option])
                else:
                    s += "{:<20} {}\n".format(option.upper(),self.opt[option])
        # atoms
        s += "\n"
        for a in self.atoms:
            s += "{:<20} {}\n".format("ATOM", a)

        # constraints
        if constr:
            for c in self.constr:
                s += "\n"
                s += "{:<20} {}\n".format("CONSTR", self._swap_constrtype(c['type']))
                for option in self.constr_keys[1:]:
                    if c[option] is not None:
                        if type(c[option])==str and "/" in c[option]:
                            s += "  {:<18} \"{}\"\n".format(option.upper(),c[option])
                        else:
                            s += "  {:<18} {}\n".format(option.upper(),c[option])
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

        # final writing
        if file is None:
            sys.stdout.write(s)
        else:
            with open(file, 'w') as f:
                f.write(s)

    def def_opt(self, mode, opt_settings):
        '''
            Assign default options for a specific mode

            Parameters
            ----------
            mode : str
                calculation mode
            opt_settings : dic
                dictionary of options in the format of 'settings.json'
        '''

        if mode is None:
            sys.stdout.write("WARNING: Default settings not assigned. Missing MODE")
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

    def resolve_constr(self, n_constr=None, step=0.05):       #TODO: take dcrd as dinit
        '''
            Resolve the iteration parameters for the
            first nth constraint: n, dinit, step & dend

            At least three are needed to resolve the fourth,
            altough 'step' can be taken from default.
            In case of conflict, 'step' is overwritten.

            Parameters
            ----------
            n_constr : int
                number of constraints to resolve, 'None' for all
            step : float
                default step lenght if not found
        '''

        # check number of constraints
        if not n_constr:
            n_constr = self.nconstr
        elif n_constr > self.nconstr:
            raise ValueError(f"More contraints requested to resolve ({n_constr}) than defined ({self.nconstr})")

        types_dict = {'dinit':float, 'dend':float, 'step':float, 'n':int}

        for nth, c in enumerate(self.constr[0:n_constr]):
            # type conversion
            for i,j in types_dict.items():
                if c[i] is not None: c[i] = j(c[i])
            # default step
            if not c['step']: c['step'] = step
            # check enough defined parameters to solve
            if (c['dinit'], c['dend'], c['n']).count(None) > 1:
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

    @staticmethod
    def _natural_sort(l):
        '''Sort a list by natural order'''
        alphanum_key = lambda key: [int(c) if c.isdigit() else c.lower() for c in re.split('([0-9]+)', key)]
        return sorted(l, key=alphanum_key)

    def launch(self):
        '''
            Launch the calculation to the queue system
        '''

        # check correct import of JACQUES
        if not jacques_import:
            sys.exit("ERROR: JACQUES could not be imported")

        # check fundamental parameters
        if self.mode is None:
            sys.exit("ERROR: Missing MODE")
        if self.opt['exe'] is None:
            sys.exit("ERROR: Missing EXE")
        if self.opt['coord'] is None:
            sys.exit("ERROR: Missing COORD")
        if self.opt['sys'] is None and self.opt['bin'] is None and self.opt['sele'] is None:
            sys.exit("ERROR: Missing topology {SYS,BIN,SELE}")

        mode   = self.mode
        name   = self.name
        out    = self.out
        opt    = self.opt
        constr = self.constr
        exe    = opt['exe']

        def submit_job(opt, jobfile):
            '''Launch to queues or display job succesfully written message'''
            if opt['dry']:
                sys.stdout.write("Job succesfully written: '{}' \n".format(jobfile))
            else:
                queues.submit(jobfile, qsys=None)

        # single structure calculations -----------------------------------
        if mode in ('sp', 'mini', 'locate', 'md', 'interaction', 'kie'):
            dynnfile = name + '.dynn'
            jobfile  = name + '.job'
            queue_param = queues.param(name, queue=opt['queue'], cores=opt['cores'], memory=opt['memory'])
            self.write(dynnfile, (mode not in ('locate')))    # no constraints for loc
            with open(jobfile, 'w') as jobf:
                jobf.write(queue_param)
                jobf.write("cd {}\n".format(os.getcwd()))
                jobf.write("{} {} > {}\n".format(exe, dynnfile, name+'.log'))
            submit_job(opt,jobfile)

        # IRC -------------------------------------------------------------
        elif mode == 'irc':
            # check direction
            if opt['irc_dir'] in (None, 0): directions = [-1,1]
            else : directions = [opt['irc_dir']]
            opt['coord'] = "../"+opt['coord']
            # loop through every direction
            for dir in directions:
                opt['irc_dir'] = dir
                # naming
                if dir == 1: name_irc = name + '-for'
                elif dir == -1: name_irc = name + '-back'
                dynnfile = name_irc + '.dynn'
                jobfile  = name_irc + '.job'
                queue_param = queues.param(name_irc, queue=opt['queue'], cores=opt['cores'], memory=opt['memory'])
                # create folder, check if hessian and move there
                shutil.rmtree(name_irc, ignore_errors=True)
                os.makedirs(name_irc)
                if os.path.isfile("update.dump"): shutil.copy("update.dump",name_irc+"/")
                os.chdir(name_irc)
                # write job, dynn and launch
                self.write(dynnfile, True)
                with open(jobfile, 'w') as jobf:
                    jobf.write(queue_param)
                    jobf.write("cd {}\n".format(os.getcwd()))
                    jobf.write("{} {} > {}\n".format(exe, dynnfile, name_irc+'.log'))
                submit_job(opt,jobfile)
                # return to workdir
                os.chdir("..")

        # POTENTIAL -------------------------------------------------------
        elif mode == 'scan':
            self.resolve_constr()
            dynnfile = name + '.dynn'
            jobfile  = name + '.job'
            queue_param = queues.param(name, queue=opt['queue'], cores=opt['cores'], memory=opt['memory'])
            self.write(dynnfile, True)
            routine = """
                      cd {pwd}\n
                      for i in $(seq 0 {n}); do
                        if [ $i == 0 ]; then
                          {exe} {dynnfile} --OUT {out} --NAME {name}.$i --N $i > {name}.$i.log
                        else
                          {exe} {dynnfile} --OUT {out} --NAME {name}.$i --N $i --COORD {name}.$il.crd > {name}.$i.log
                        fi
                        il=$i
                      done
                      """.format(pwd=os.getcwd(), n=constr[0]['n'], exe=exe, dynnfile=dynnfile, name=name, out=out)
            with open(jobfile, 'w') as jobf:
                jobf.write(queue_param)
                jobf.write(dedent(routine))
            submit_job(opt,jobfile)

        elif mode == 'pes':
            self.resolve_constr()
            dynnfile = name + '.dynn'
            self.write(dynnfile, True)
            # first constraint preference
            for i in range(0, int(constr[0]['n'])+1):
                name_pes = "{}.{}".format(name, i)
                jobfile  = name_pes + '.job'
                jobn     = "{}.{}.job".format(name, i+1)
                queue_param = queues.param(name_pes, queue=opt['queue'], cores=opt['cores'], memory=opt['memory'])
                if i==0: coord0 = opt['coord']
                else: coord0 = "{}.{}.$j.crd".format(name, i-1)
                routine = """
                          cd {pwd}\n
                          for j in $(seq 0 {n}); do
                            if [ $j == 0 ]; then
                              {exe} {dynnfile} --OUT {out}.{i}.out --NAME {name}.{i}.$j --N {i} $j --COORD {coord0} > {name}.{i}.$j.log
                              {{ qsub   {jobn} ; }} 2>/dev/null
                              {{ sbatch {jobn} ; }} 2>/dev/null
                            else
                              {exe} {dynnfile} --OUT {out}.{i}.out --NAME {name}.{i}.$j --N {i} $j --COORD {name}.{i}.$jl.crd > {name}.{i}.$j.log
                            fi
                            jl=$j
                          done
                          """.format(pwd=os.getcwd(), n=constr[1]['n'], exe=exe,
                                     dynnfile=dynnfile, name=name, out=out.split(".out")[0], i=i, jobn=jobn, coord0=coord0)
                with open(jobfile, 'w') as jobf:
                    jobf.write(queue_param)
                    jobf.write(dedent(routine))
            submit_job(opt, name+".0.job")

        # FREE ENERGY / CORRECTION ------------------------------------
        elif mode in ('pmf', 'corr'):
            self.resolve_constr(self.dim)
            dynnfile = name + '.dynn'
            self.write(dynnfile, True)
            # create jobs folder
            shutil.rmtree("jobs", ignore_errors=True)
            os.mkdir("jobs")
            # get list of crd files
            crd_dir = opt['coord']
            if not os.path.isdir(crd_dir):
                sys.exit("ERROR: Path to look for crd not found")
            crd_files = self._natural_sort(glob.glob(os.path.join(crd_dir, "*.crd")))
            # loop through all requested dimensions
            job_files = []
            not_found = []
            for c in itertools.product(*[range(constr[i]['n']+1) for i in range(self.dim)]):
                numd = ".".join(map(str, c))
                nums = " ".join(map(str, c))
                jobfile = name + "." + numd + ".job"
                queue_param = queues.param(name+"."+numd, queue=opt['queue'], cores=opt['cores'], memory=opt['memory'])
                # find corresponding coordinates
                try:
                    crd = [i for i in crd_files if "."+numd+".crd" in i][0]
                except IndexError:
                    not_found.append(nums)
                    continue
                if mode == 'pmf':
                    routine = """
                              cd {pwd}\n
                              {exe} {dynnfile} --NAME {name}.{numd} --N {nums} --COORD {crd} > {name}.{numd}.log\n
                              """.format(pwd=os.getcwd(), exe=exe, dynnfile=dynnfile, name=name, numd=numd, nums=nums, crd=crd)
                elif mode == 'corr':
                    routine = """
                              mkdir {pwd}/{name}.{numd}
                              cd {pwd}/{name}.{numd}\n
                              {exe} ../{dynnfile} --OUT ../{out} --NAME {name}.{numd} --N {nums} --COORD ../{crd} > {name}.{numd}.log\n
                              """.format(pwd=os.getcwd(), exe=exe, dynnfile=dynnfile, name=name, out=out, numd=numd, nums=nums, crd=crd)
                with open(os.path.join("jobs", jobfile), 'w') as jobf:
                    jobf.write(queue_param)
                    jobf.write(dedent(routine))
                job_files.append(jobfile)
            # not found warning
            if not_found:
                with open("crd_not_found.txt", 'w') as f:
                    f.write("\n".join(not_found))
                sys.stdout.write("WARNING: Some crd could not be found. Registered on 'crd_not_found.txt'.\n")
            # main job launcher driver
            with open(name+".job", 'w') as jobf:
                jobf.write(queues.param(name, queue=opt['queue'], cores=1))
                jobf.write(f"\ncd {os.getcwd()}\n\n")
                jobf.write("jobs=(\\\n"+"\n".join(job_files)+"\n)\n")
                routine = r"""
                          for job in ${jobs[@]}; do
                              { qsub   jobs/$job ; } 2>/dev/null
                              { sbatch jobs/$job ; } 2>/dev/null
                              sleep 1
                          done
                          """
                jobf.write(dedent(routine))
            submit_job(opt, name+".job")

        # UNKNOWN ---------------------------------------------------------
        else:
            sys.exit(f"ERROR: Unkown mode '{mode}'")

    def post(self, rm=True):
        '''
            Post-process routine after a calculation

            Parameters
            ----------
            rm : bool
                remove files after processing
        '''

        filetypes = ('log', 'crd', 'out', 'dat')

        mode   = self.mode
        name   = self.name

        # check fundamental parameters
        if mode is None:
            sys.exit("ERROR: Missing MODE")
        elif mode not in ('scan', 'pes', 'pmf', 'corr'):
            sys.exit(f"ERROR: No post-process routine for this mode '{mode}'")

        # list of files to process in natural order
        dir_files = dict()
        for filetype in filetypes:
            if filetype in ('log', 'crd', 'out', 'job'):
                file_pattern = f"*.{filetype}"
            elif filetype == 'dat':
                file_pattern = f"{filetype}_*"
            dir_files[filetype] = self._natural_sort(glob.glob(file_pattern))

        sys.stdout.write(f"## POST-PROCESS: {mode}\n# NAME: {name}\n\n")

        for filetype in filetypes:
            # check no files matching
            if dir_files[filetype]:
                n_tot = len(dir_files[filetype])
                # pre-loop actions
                if filetype in ('log', 'out'):
                    final_file = f"{name}.{filetype}"
                    # remove final file from file lists
                    if final_file in dir_files[filetype]: dir_files[filetype].remove(final_file)
                    f_final = open(final_file, 'a+')
                elif filetype == 'crd' and not os.path.isdir("crd"):
                    os.mkdir("crd")
                elif filetype == 'dat' and not os.path.isdir("dat"):
                    os.mkdir("dat")
                # loop through every matched file
                for n, f in enumerate(dir_files[filetype]):
                    # LOG / OUT -----
                    if filetype in ('log', 'out'):
                        if filetype == 'log': f_final.write(f"\n#######  {f}  \n")
                        with open(f, 'r') as f_tmp:
                            f_final.write(f_tmp.read())
                        if rm: os.remove(f)
                    # CRD -----
                    elif filetype == 'crd':
                        shutil.move(f, os.path.join("crd", f))
                    # DAT -----
                    elif filetype == 'dat':
                        shutil.move(f, os.path.join("dat", f))
                    # progress bar
                    sys.stdout.write("\r# {}s  [{:21s}] - {:>6.2f}%".format(filetype.upper(),
                                                                            "■"*int(21.*(n+1)/n_tot),
                                                                            (n+1)/n_tot*100))
                    sys.stdout.flush()
                sys.stdout.write(f"\r# {filetype.upper()}s  ✔ {' ':<40}\n")
                # post-loop actions
                if filetype in ('log', 'out'):
                    f_final.close()
            else:
                sys.stdout.write(f"# {filetype.upper()}s  ⨯\n")

        # JOB -----
        if rm and os.path.isdir("jobs"):
            shutil.rmtree("jobs")
            sys.stdout.write("# JOBSs  ✔\n")

        sys.stdout.write("\n")
