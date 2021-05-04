"""
=======================================================================
  Main execution routine
=======================================================================

  Functions
  ---------

    main

"""

import sys

from jacques                import queues_def, opt_def, dynamon_path
from jacques._parser        import parser
from jacques.dynnconfig     import DynnConfig


def main():
    # check no argument
    if len(sys.argv) <= 1:
        sys.exit("ERROR: No input arguments. Use -h for help.")

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
