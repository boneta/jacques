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
from jacques.general_parser import parser
from jacques.dynnconfig     import DynnConfig


def main():
    # check no argument
    if len(sys.argv) <= 1:
        sys.exit("ERROR: No input arguments. Use -h for help.\n")

    # parser
    p = parser(sys.argv[1])
    args = p.parse_args()

    # configuration object
    config = DynnConfig()
    if args.f is not None:
        config.read_file(args.f)
    for option in vars(args):
        value = getattr(args, option)
        if value == False: value = None
        if option in DynnConfig.opt_keys and value is not None:
            config.opt[option] = value
    config.def_opt(config.mode, opt_def)

    # configure and launch calculation
    config.launch()

if __name__ == '__main__':
    main()
