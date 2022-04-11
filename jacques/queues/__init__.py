"""
=======================================================================
  Queues and job management
=======================================================================

Queues
------

    QueueManager

Job files
---------

    JobFile

"""

from jacques import jobf, msgf, settings


default_queues = settings['queues']
supported_managers = {'slurm': 'sbatch',
                        'sge': 'qsub'}

from .jobfile import *
from .manager import *
