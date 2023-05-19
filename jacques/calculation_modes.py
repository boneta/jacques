"""
=======================================================================
  Calculation routines
=======================================================================

  Functions
  ---------

    launch_calc

"""

import glob
import os
import re
import sys
from copy import deepcopy
from textwrap import dedent

from .dynnconfig import DynnConfig
from .queues import JobFile


def _natural_sort(l:list) -> list:
    '''Sort a list by natural order'''
    alphanum_key = lambda key: [int(c) if c.isdigit() else c.lower() for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)

def launch_calc(dynnconfig:'DynnConfig') -> None:
    '''
        Launch the calculation to the queue system

        Parameters
        ----------
        dynnconfig: DynnConfig
            DynnConfig object
    '''

    config = deepcopy(dynnconfig)

    # check fundamental parameters
    if config.mode is None:
        sys.exit("ERROR: Missing MODE")
    if config.opt['exe'] is None:
        sys.exit("ERROR: Missing EXE")
    if config.opt['coord'] is None:
        sys.exit("ERROR: Missing COORD")
    if config.opt['sys'] is None and config.opt['bin'] is None and config.opt['sele'] is None:
        sys.exit("ERROR: Missing topology {SYS,BIN,SELE}")

    match config.mode:

        # SINGLE STRUCTURE --------------------------------------------
        case 'sp' | 'mini' | 'md' | 'interaction' | 'kie' | 'locate':
            if config.mode == 'locate':    # no constraints for locate
                config.constr = []
            config.write_dynn()
            routine = f"{config.opt['exe']} {config.dynn} > {config.name}.log\n"

        # IRC ---------------------------------------------------------
        case 'irc':
            # both directions (0)
            if config.opt['irc_dir'] in (None, 0, '0'):
                config.irc_both_dir = True
                config.opt['coord'] = f"../{config.opt['coord']}"
                config.opt['out'] = f"../{config.out}"
                config.write_dynn()
                name = config.name
                config.dynn0 = config.dynn
                for dir in [-1, 1]:
                    config.opt['irc_dir'] = dir
                    config.opt['name'] = f"{name}-BACK" if dir == -1 else f"{name}-FOR"
                    launch_calc(config)
                return
            # single direction (-1/1)
            else:
                # comes from both direction call
                # (create folder, move there, copy hessian if present, run)
                if hasattr(config, 'irc_both_dir'):
                    routine = f"mkdir -p {config.name}\ncd {config.name}\n"
                    if os.path.isfile("update.dump"):
                        routine += f"cp ../update.dump .\n"
                    routine += f"{config.opt['exe']} ../{config.dynn0} --NAME {config.name} --IRC_DIR {config.opt['irc_dir']} > {config.name}.log\n"
                # normal (just run)
                else:
                    config.write_dynn()
                    routine = f"{config.opt['exe']} {config.dynn} > {config.name}.log\n"

        # POTENTIAL ---------------------------------------------------
        case 'scan':
            config.resolve_constr(1)
            config.write_dynn()
            routine = """
                    set -e
                    coord={coord0}
                    for i in {{0..{n}}}; do
                        {exe} {dynnfile} --NAME {name}.$i --N $i --COORD $coord >> {name}.log
                        coord={name}.$i.crd
                    done
                    """.format(name=config.name,
                               dynnfile=config.dynn,
                               exe=config.opt['exe'],
                               coord0=config.opt['coord'],
                               n=config.constr[0]['n'])

        case 'pes':
            config.resolve_constr(2)
            config.write_dynn()
            config.opt['array_first'] = 0
            config.opt['array_last'] = config.constr[0]['n']
            routine = """
                    set -e
                    if [ $ID -eq 0 ]; then
                        coord={coord0}
                    else
                        coord={name}.$(($ID-1)).0.crd
                    fi
                    for waiting in {{1..360}}; do
                        [[ -f $coord ]] && break
                        sleep 10s
                    done
                    for j in {{0..{n}}}; do
                        {exe} {dynnfile} --NAME {name}.$ID.$j --N $ID $j --COORD $coord --OUT {name}.$ID.out >> {name}.$ID.log
                        coord={name}.$ID.$j.crd
                    done
                    """.format(name=config.name,
                               dynnfile=config.dynn,
                               exe=config.opt['exe'],
                               coord0=config.opt['coord'],
                               n=config.constr[1]['n'])

        # FREE ENERGY / CORRECTION ------------------------------------
        case 'pmf' | 'corr':
            config.resolve_constr(config.dim)
            # get list of crd files
            crd_dir = config.opt['coord']
            if not os.path.isdir(crd_dir):
                sys.exit("ERROR: Path to look for crd not found")
            crd_files = _natural_sort(
                glob.glob(os.path.join(crd_dir, "*.crd")))
            crd_files = [os.path.basename(crd) for crd in crd_files]
            # correction path
            dynnfile = config.dynn
            if config.mode == 'corr':
                dynnfile = f"../{dynnfile}"
                config.opt['coord'] = f"../{config.opt['coord']}"
                config.opt['out'] = f"../{config.out}"
            config.write_dynn()
            # build argument for each crd file (based on dot separated numbers in crd files)
            arguments = [
                f"'{config.name}  {' '.join(crd_files[-1].split('.')[1:-1])}'"]
            for crd in crd_files:
                num = crd.split('.')[1:-1]
                arguments.append("'--NAME {name}.{num_dot} --N {num_spc} --COORD {crd_path}'" \
                                .format(name=config.name,
                                        num_dot='.'.join(num),
                                        num_spc=' '.join(num),
                                        crd_path=os.path.join(config.opt['coord'], crd)))
            # build routine
            routine = "\narg=( {} )\n\n".format(' \\\n      '.join(arguments))
            if config.mode == 'corr':
                routine += f"mkdir -p {config.name}.$ID\ncd {config.name}.$ID\n"
            routine += "{exe} {dynnfile} ${{arg[$ID]}} > {name}.${{ID}}.log\n" \
                    .format(name=config.name,
                            dynnfile=dynnfile,
                            exe=config.opt['exe'])
            # set-up array
            config.opt['array_first'] = 1
            config.opt['array_last'] = len(crd_files)

        # UNKNOWN -----------------------------------------------------
        case _:
            sys.exit(f"ERROR: Unkown mode '{config.mode}'")

    jobfile = JobFile(dedent(routine), **config.opt)
    jobfile.submit(config.opt['dry'])
