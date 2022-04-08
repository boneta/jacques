"""
=======================================================================
  Cluster queues management
=======================================================================

  Functions
  ---------

    param
    submit
    _conv_mem

"""

import re
import subprocess
import sys
from textwrap import dedent

from jacques import jobf, msgf, queues_def

from . import submit_exe, queue_templates


def param(name, queue=None, cores=1, memory='4000MB', memory_plus='1000MB',
          nodes=1, queues_def=queues_def, msgf=msgf, jobf=None):
    '''
        Build queue parameters in bash format

        Parameters
        ----------
        name : str
            name for the job
        queue : {None, 'sge', 'slurm'}, optional
            queue name. if None, used from queues_def
        cores : int, optional
            number of total CPU cores to request (per node)
        memory : str, optional
            amount of total RAM memory to request (per node)
        memory_plus : str, optional
            amount of extra RAM memory to request
        nodes : int, optional
            number of nodes to request
        queues_def : dic
            dictionary of default queues {'queue_sys':'queue_name'}
        msgf : str, optional
            path to redirect output and error messages from queue system
        jobf : str, optional
            path to move the script file after execution, skipped if None

        Returns
        -------
        queue_param : str
            string of queue parameters in bash format

    '''

    if queue is None:
        queue_sge = queues_def['sge']
        queue_slurm = queues_def['slurm']
    else:
        queue_sge = queue
        queue_slurm = queue

    if memory is None: memory = '4000MB'
    memory_final = "{}MB".format(round(_conv_mem(memory,'MB')+_conv_mem(memory_plus,'MB')))

    # Sun Grid Engine
    sge = queue_templates['SGE'].format(name=name,
                                        queue=queue_sge,
                                        msgf=msgf,
                                        cores=cores)
    if cores == 1:
        sge = sge.replace('#$ -pe mp1 1\n', '')

    # SLURM
    slurm = queue_templates['SLURM'].format(name=name,
                                            queue=queue_slurm,
                                            msgf=msgf,
                                            nodes=nodes,
                                            cores=cores,
                                            memory=memory_final)

    # move job file
    mvjob = "\nmv $0 {}\n".format(jobf) if jobf is not None else "\n"

    queue_param = "#!/bin/bash\n\n" + dedent(sge) + "\n"+ dedent(slurm) + "\n" + mvjob

    return queue_param

def submit(jobfile, submit_exe=submit_exe):
    '''
        Submit a job file to the queue system

        Parameters
        ----------
        jobfile : str
            name of the job file
        submit_exe : str, optional
            name of the submitter executable (i.e. sbatch/qsub)
            if None, a message is shown and not submitted
            if omitted is guessed from the environment
    '''

    if submit_exe is not None:
        subprocess.call([submit_exe, jobfile])
    else:
        sys.stdout.write(f"ERROR: No queue system detected. Job '{jobfile}' not submitted\n")

def _conv_mem(memory, out_units='MB'):
    '''
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
    '''

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
