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

from jacques import dynamon_path, jacques_path, opt_def
from jacques._parser import parser
from jacques.dynnconfig import DynnConfig


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
      config.post()
    else:
      config.launch()

if __name__ == '__main__':
    main()
