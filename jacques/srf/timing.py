"""
=======================================================================
  Timing
=======================================================================

Classes
-------

    ElapTime

"""

import time

class ElapTime:
    """
        Elapsed time accounting class

        Attributes
        ----------
        t0 : float

        Properties
        ----------
        elapsed : float
            elapsed time since start (s)

        Methods
        -------
        reset()
            reset the timer and start counting from now

        Static Methods
        --------------
        now_fmt() -> str
            get the current date-time in complete format
            e.g.: 'Wed Jun  9 04:26:40 1993'

    """

    def __init__(self):
        """ Start counting time since epoch """
        self.t0 = time.time()

    def __repr__(self):
        dt = self.elapsed
        return f"<ElapTime with {dt}>"

    def __str__(self):
        dt = self.elapsed
        return f"{int(dt//86400)}-{time.strftime('%H:%M:%S', time.gmtime(dt))}"

    def reset(self):
        """ Reset the timer and start counting from now """
        self.__init__()

    @property
    def elapsed(self) -> float:
        """ Elapsed time since start """
        return time.time() - self.t0

    @staticmethod
    def now_fmt() -> str:
        """ Formatted current date-time """
        return time.ctime(time.time())
