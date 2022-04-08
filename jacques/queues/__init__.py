"""
=======================================================================
  Queue management
=======================================================================

Queues
------

    param
    submit

"""

import os
from shutil import which


# read queue templates
queue_templates = {'SGE': "", 'SLURM': ""}
for queue_sys in queue_templates.keys():
    with open(os.path.join(os.path.dirname(__file__), "templates", f"header_{queue_sys}.sh"), 'r') as f:
        queue_templates[queue_sys] = f.read()

# detect queue system and submit executable
submit_exe = None
for exe in ['sbatch', 'qsub']:
    if submit_exe := which(exe):
        break


from .queues import *
