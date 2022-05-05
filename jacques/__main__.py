"""
=======================================================================
  Main execution routine
=======================================================================

  Functions
  ---------

    main

"""

import os
import subprocess
import sys

from . import dynamon_path, jacques_path, opt_def
from ._parser import parser
from .dynnconfig import DynnConfig
from .processing import post


def main():
    # no arguments
    if len(sys.argv) <= 1:
        sys.argv.append('-h')
    # special arguments
    if sys.argv[1] == 'update':
        subprocess.run([os.path.join(jacques_path, "update.sh")])
        sys.exit(0)

    # parser
    args = parser(sys.argv[1])

    # configuration object
    config = DynnConfig()
    config.read_opt(**vars(args))
    config.def_opt(config.mode, opt_def)

    # run post-process or launch calculation
    if config.opt['post']:
      post(config)
    else:
      config.launch()

if __name__ == '__main__':
    main()
