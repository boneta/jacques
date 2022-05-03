"""
=======================================================================
  Cluster queues management
=======================================================================

  Classes
  -------

    JobFile

"""

import re
import sys
from copy import deepcopy

from . import jobf, msgf
from .manager import QueueManager


class JobFile:
    """
        Job file class

        Includes resources to be requested and queue parameters;
        and subroutines/instructions to do the job.

        Parameters
        ----------
        name : str, optional
            name for the job (def: 'job')
        commands : str, optional
            formatted sring of commands to be executed as job
        queue_manager : str, optional
            queue manager to force use (e.g. slurm/sge)
            by default it will be guessed from the system PATH
        kwargs : dict, optional
            additional queue parameters

        Attributes
        ----------
        name : str
            name for the job (def: 'job')
        jobfile : str
            name of the job file (def: 'job.job')
        commands : str
            formatted sring of commands to be executed as job
        opt : dict
            additional queue parameters
            name, queue, cores, memory, memory_plus, nodes, array_first, array_last, msgf, jobf

        Methods
        -------
        write()
            Write to a file (.job)
        submit(dry)
            Write and submit the job file to the queue manager
    """

    opt_def = {
        'name': 'job',
        'queue': '',
        'cores': 1,
        'memory': '4000MB',
        'memory_plus': '1000MB',
        'nodes': 1,
        'array_first': 0,
        'array_last': 0,
        'msgf': msgf,
        'jobf': None
        }

    def __init__(self, commands="", queue_manager=None, **kwargs) -> None:
        """
            Build queue parameters

            Parameters
            ----------
            commands : str, optional
                formatted sring of commands to be executed as job
            queue_manager : str, optional
                queue manager to force use (e.g. slurm/sge)
                by default it will be guessed from the system PATH
            kwargs : dict, optional
                additional queue parameters:
                name, queue, cores, memory, memory_plus, nodes, array_first, array_last, msgf, jobf
        """
        self.commands = commands
        self.qmanager = QueueManager(queue_manager)
        self.opt = deepcopy(self.opt_def)
        self.opt['queue'] = self.qmanager.queue_def
        for key, value in self.opt.items():
            self.opt[key] = kwargs.get(key, value) or value

    @property
    def name(self) -> str:
        return self.opt['name']

    @property
    def jobfile(self) -> str:
        return f"{self.name}.job"

    def __str__(self) -> str:
        """
            Build queue parameters as a bash script

            Returns
            -------
            queue_param : str
                string of queue parameters in bash format
        """
        # modify some parameters
        opt = deepcopy(self.opt)
        final_memory = round(_conv_mem(opt['memory'], 'MB') + _conv_mem(opt['memory_plus'], 'MB'))
        opt['memory'] = f"{final_memory}MB"
        # queue manager options and resources
        queue_text = ""
        if self.qmanager:
            queue_text = self.qmanager.template.format(**opt)
            # remove parallel environment for one core in SGE and array if only one job
            queue_text = queue_text.replace('#$ -pe mp1 1\n', '') \
                                   .replace('#$ -t 0-0\n', '') \
                                   .replace('#SBATCH --array=0-0\n', '')
        # combined text
        return f"#!/bin/bash\n\n{queue_text}\n{self.commands}\n"

    def write(self) -> None:
        """Write to a file (.job)"""
        with open(self.jobfile, 'w') as f:
            f.write(str(self))

    def submit(self, dry=False) -> None:
        """
            Write and submit the job file to the queue manager

            Parameters
            ----------
            dry : bool, optional
                dry run (def: False)
                display a message but do not submit the job
        """
        self.write()
        if not dry:
            self.qmanager.submit(self.jobfile)
            sys.stdout.write(f"Jobfile submitted: {self.jobfile}\n")
        else:
            sys.stdout.write(f"Jobfile written: {self.jobfile}\n")


def _conv_mem(memory, out_units='MB'):
    """
        Convert between RAM memory formats

        Parameters
        ----------
        memory : str / float
            memory to be converted
            string format (i.e. '2000MB') or float (<100: GB, >100:MB)
        out_units : {KB, MB, GB, TB, MW, GW}, optional
            memory units for output

        Returns
        -------
        memory_conv : float
            memory converted in out_units
    """

    # memory units with MB as base
    units = {'KB':1.0E-3, 'MB':1.0, 'GB':1.0E3, 'TB':1.0E6, 'MW':8.0, 'GW':8.0E3}

    # assign amount of memory and input units
    try:
        mem = float(re.findall(r'[0-9]*\.?[0-9]+', memory)[0])
        in_units = re.findall(r'[A-Za-z]+', memory)[0].upper()
    except:
        mem = float(memory)
        in_units = 'GB' if mem < 100 else 'MB'  # guess units (risky)

    out_units = out_units.upper()

    if in_units  not in units.keys(): raise NameError('Wrong input memory units')
    if out_units not in units.keys(): raise NameError('Wrong output memory units')

    # convert
    memory_conv = mem * units[in_units] / units[out_units]

    return memory_conv
