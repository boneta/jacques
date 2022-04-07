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


# read queue templates
queue_templates = {'SGE': "", 'SLURM': ""}
for queue_sys in queue_templates.keys():
    with open(os.path.join(os.path.dirname(__file__), "templates", f"header_{queue_sys}.sh"), 'r') as f:
        queue_templates[queue_sys] = f.read()


from .queues import *
