"""
=======================================================================
  Main execution routine
=======================================================================

  Functions
  ---------

    main

"""

import os
import sys
import subprocess

from jacques                import queues_def, opt_def, jacques_path, dynamon_path
from jacques._parser        import parser
from jacques.dynnconfig     import DynnConfig


def main():
    # special arguments -> none or update
    if len(sys.argv) <= 1:
        sys.exit("ERROR: No input arguments. Use -h for help.")
    elif sys.argv[1] == 'update':
        subprocess.run([os.path.join(jacques_path, "update.sh")])
        sys.exit(0)

    # parser
    p = parser(sys.argv[1])
    args = p.parse_args()

    # configuration object
    config = DynnConfig()
    config.read_opt(**vars(args))
    config.def_opt(config.mode, opt_def)

    # configure and launch calculation
    config.launch()

if __name__ == '__main__':
    main()
