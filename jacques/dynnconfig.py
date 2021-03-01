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
import shutil
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
        opt : list
        atoms : list of list
        constr : list of dict

        Properties
        ----------
        mode : str
            calculation mode
        name : str
            calculation name
        nconstr : int
            number of defined constraints

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
    '''

    opt_keys = [
        'mode',

        'name', 'sys', 'bin', 'sele', 'coord', 'ff', 'seq',
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

        'exe', 'queue', 'jobonly'
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
        if file is not None: self.read_file(file)

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
    def nconstr(self):
        '''Number of defined constraints'''
        return len(self.constr)

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
            opts = list(map(str.strip, opts))  # remove trailing spaces

        # remove comment lines or strange empty lines and space-split the first word
        opts = [list(map(str.strip,line.split(' ',1))) for line in opts
                if len(line) > 0 and not line.startswith("#") and not line.startswith("!")]

        # loop line by line to assign options
        n = 0
        while n < len(opts):
            line = list(map(str.strip,opts[n]))
            if line[0] in (self.opt_keys_up|self.opt_keys_low):
                self.opt[line[0].lower()] = line[1].strip('"')
            elif line[0] == 'ATOM' or line[0] == 'A':
                self.atoms.append(line[1])
            elif line[0] == 'CONSTR' or line[0] == 'C':
                self.constr.append(self.constr_empty.copy())
                self.constr[-1]['type'] = self._swap_constrtype(line[1])      # type of constr
                # read constraint options
                while n < len(opts):
                    n += 1
                    line = opts[n]
                    if line[0] == 'CONSTR' or line[0] == 'C': break
                    if line[0] in (self.constr_keys_up|self.constr_keys_low):
                        self.constr[-1][line[0].lower()] = line[1].strip('"')
                    else:
                        print("WARNING: Unrecognized option:", line[0])
            else:
                print("WARNING: Unrecognized option:", line[0])
            n += 1

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

        if file is None:
            f = sys.stdout
        else:
            f = open(file, 'wt')

        # mode
        try:
            f.write("{:<20} {}\n\n".format("MODE",self.opt['mode'].upper()))
        except AttributeError:
            sys.stdout.write("WARNING: 'MODE' parameter not specified\n")

        # general options
        for option in self.opt_keys[1:-3]:
            if self.opt[option] is not None:
                if type(self.opt[option])==str and ("/" in self.opt[option] or "," in self.opt[option]):
                    f.write("{:<20} \"{}\"\n".format(option.upper(),self.opt[option]))
                else:
                    f.write("{:<20} {}\n".format(option.upper(),self.opt[option]))
        # atoms
        f.write("\n")
        for a in self.atoms:
            f.write("{:<20} {}\n".format("ATOM", a))

        # constraints
        if constr:
            for c in self.constr:
                f.write("\n")
                f.write("{:<20} {}\n".format("CONSTR", self._swap_constrtype(c['type'])))
                for option in self.constr_keys[1:]:
                    if c[option] is not None:
                        if type(c[option])==str and "/" in c[option]:
                            f.write("  {:<18} \"{}\"\n".format(option.upper(),c[option]))
                        else:
                            f.write("  {:<18} {}\n".format(option.upper(),c[option]))
                f.write("C\n")

        if file is not None: f.close()

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

        for nth, c in enumerate(self.constr[0:n_constr]):
            if not c['step']: c['step'] = step
            # check enough defined parameters to solve
            if (c['dinit'], c['dend'], c['n']).count(None) > 1:
                raise ValueError(f"Not enough parameters specified to resolve constraint {nth+1}")
            # dinit & dend -> n [& step]
            if c['dinit'] and c['dend']:
                diff = float(c['dend']) - float(c['dinit'])
                if not c['n']:
                    c['n'] = m.ceil(diff/c['step'])
                    if c['n'] < 0: c['n'], c['step'] = -c['n'], -c['step']
                else:
                    c['step'] = round(diff/c['n'], 3)
            # dinit & n [& step] -> dend
            elif c['dinit']:
                c['dend'] = float(c['dinit']) + float(c['step'])*int(c['n'])
            # dend & n [& step] -> dinit
            elif c['dend']:
                c['dinit'] = float(c['dend']) - float(c['step'])*int(c['n'])

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
        opt    = self.opt
        constr = self.constr
        exe    = opt['exe']

        def submit_job(opt, jobfile):
            '''Launch to queues or display job succesfully written message'''
            if opt['jobonly']:
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
        elif mode in ('irc'):
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
        elif mode in ('scan'):
            self.resolve_constr()
            dynnfile = name + '.dynn'
            jobfile  = name + '.job'
            queue_param = queues.param(name, queue=opt['queue'], cores=opt['cores'], memory=opt['memory'])
            self.write(dynnfile, True)
            routine = """
                      cd {pwd}\n
                      for i in $(seq 0 {n}); do
                        if [ $i == 0 ]; then
                          {exe} {dynnfile} --NAME scan.$i --N $i > scan.$i.log
                        else
                          {exe} {dynnfile} --NAME scan.$i --N $i --COORD scan.$il.crd > scan.$i.log
                        fi
                        il=$i
                      done
                      """.format(pwd=os.getcwd(), n=constr[0]['n'], exe=exe, dynnfile=dynnfile)
            with open(jobfile, 'w') as jobf:
                jobf.write(queue_param)
                jobf.write(dedent(routine))
            submit_job(opt,jobfile)

        elif mode in ('pes'):
            self.resolve_constr()
            dynnfile = name + '.dynn'
            self.write(dynnfile, True)
            # first constraint preference
            for i in range(0, int(constr[0]['n'])+1):
                name_pes = "pes.{}".format(i)
                jobfile  = name_pes + '.job'
                jobn     = "pes.{}.job".format(i+1)
                queue_param = queues.param(name_pes, queue=opt['queue'], cores=opt['cores'], memory=opt['memory'])
                if i==0: coord0 = opt['coord']
                else: coord0 = "pes.{}.$j.crd".format(i-1)
                routine = """
                          cd {pwd}\n
                          for j in $(seq 0 {n}); do
                            if [ $j == 0 ]; then
                              {exe} {dynnfile} --NAME pes.{i} --N {i} $j --COORD {coord0} > pes.{i}.$j.log
                              {{ qsub   {jobn} ; }} 2>/dev/null
                              {{ sbatch {jobn} ; }} 2>/dev/null
                            else
                              {exe} {dynnfile} --NAME pes.{i} --N {i} $j --COORD pes.{i}.$jl.crd > pes.{i}.$j.log
                            fi
                            jl=$j
                          done
                          """.format(pwd=os.getcwd(), n=constr[1]['n'], exe=exe,
                                     dynnfile=dynnfile, i=i, jobn=jobn, coord0=coord0)
                with open(jobfile, 'w') as jobf:
                    jobf.write(queue_param)
                    jobf.write(dedent(routine))
            submit_job(opt,"pes.0.job")

        # FREE ENERGY -----------------------------------------------------
        elif mode in ('pmf'):
            sys.exit("Mode not implemented yet")
        # CORRECTION ------------------------------------------------------
        elif mode in ('corr'):
            sys.exit("Mode not implemented yet")
        # UNKNOWN ---------------------------------------------------------
        else:
            sys.exit(f"ERROR: Unkown mode '{mode}'")
