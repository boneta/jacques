"""
=======================================================================
  Queue manager
=======================================================================

Classes
-------

    QueueManager

"""

import os.path
import shlex
import sys
from shutil import which
from subprocess import call

from . import default_queues, supported_managers


class QueueManager:
    """
        Queue manager class

        Parameters
        ----------
        manager : str, optional
            queue manager (e.g. slurm/sge)

        Attributes
        ----------
        manager : str
            queue manager (e.g. slurm/sge)
        exe : str
            submit executable (e.g. sbatch/qsub)
        queue_def : str
            default queue/partition name
        template : str
            template for the job file

        Methods
        -------
        guess_manager()
            Try to detect the queue manager available
        submit(jobfile)
            Submit a job file to the queue manager
    """

    supported_managers = supported_managers

    def __init__(self, manager=None) -> None:
        """
            Initialize the queue manager

            A user defined queue manager can be specified,
            otherwise it will be guessed from the system PATH.

            Assign submit executable, read the queue template
            and default queue/partition name based on the queue manager

            Parameters
            ----------
            manager : str, optional
                queue manager (e.g. slurm/sge)
        """
        if manager is not None and manager not in supported_managers:
            raise ValueError(f"Unknown queue manager: {manager}")
        self.manager = manager or self.guess_manager()
        if self.manager:
            self.exe = self.supported_managers[self.manager]
            self.queue_def = default_queues[self.manager]
            with open(os.path.join(os.path.dirname(__file__), "templates", f"header_{self.manager.upper()}.sh"), 'r') as f:
                self.template = f.read()
        else:
            self.exe, self.queue_def, self.template = None, "", None

    def __bool__(self) -> bool:
        return bool(self.manager)

    @staticmethod
    def guess_manager() -> str:
        """
            Try to detect the queue manager available

            Returns
            -------
            manager : str/None
                queue manager (e.g. slurm/sge)
                None if no queue manager is found
        """
        for manager, exe in supported_managers.items():
            if which(exe):
                return manager
        else:
            return None

    def submit(self, jobfile, options="", dry=False) -> None:
        """
            Submit a job file to the queue manager

            Parameters
            ----------
            jobfile : str
                job file (bash script)
            options : str, optional
                additional options to pass to the queue manager
            dry : bool, optional
                dry run (do not submit the job)
        """
        submit_command = [self.exe, jobfile]
        if options:
            submit_command = [submit_command[0], *shlex.split(options, posix=False), *submit_command[1:]]
        if self:
            if dry:
                sys.stdout.write(shlex.join(submit_command) + "\n")
            else:
                call(submit_command)
        else:
            sys.stdout.write(f"ERROR: No queue system detected. Job '{jobfile}' not submitted\n")
