"""
=======================================================================
  Parser
=======================================================================

  Functions
  ---------

    parser

"""

import argparse

from jacques import __version__


def parser(mode=None):
    '''
        Construct specific parser for a MODE

        Parameters
        ----------
        mode : str
            calculation mode

        Returns
        -------
        p : ArgumentParser
            argparse's parser object
    '''

    if mode: mode = mode.lower()
    mode_list = ('sp', 'mini', 'locate', 'md', 'interaction',
                 'kie', 'irc', 'pes', 'scan', 'pmf', 'corr')

    # custom argparse to include help
    class MyParser(argparse.ArgumentParser):
        def format_help(self):
             return h

    # help text
    h=  " ##################  JACQUES  ##################\n"
    h+=" -- Friendly interface for QM/MM calculations --\n\n"

    if mode not in mode_list: mode = "<mode>"
    h+=f" USAGE:   jacques {mode} [[--option arg] ...]\n\n"
    h+=" OPTIONS:\n"

    p = MyParser(formatter_class=argparse.RawTextHelpFormatter, add_help=False)

    #  MODE  ----------------------------------------------------------
    if mode not in mode_list:
        h+=""" <mode>                            main calculation mode
                                     sp            single point
                                     mini          minimization/optimization
                                     locate        locate and characterize minimum/TS
                                     md            molecular dynamics
                                     interaction   electrostatic/VdW interactions along trajectory
                                     kie           kinetic isotope effects
                                     irc           internal reaction coordinate from TS
                                     scan          potential energy surface (1D)
                                     pes           potential energy surface (2D)
                                     pmf           potential of the mean force
                                     corr          spline corrections
            \n"""
            
    p.add_argument('mode', type=str, nargs='?', choices=mode_list)

    #  CONFIG FILE  ---------------------------------------------------
    h+=" -f  <.jcq>                        file config in DYNAMON format\n\n"
    p.add_argument('-f', type=str, required=(mode in ('scan', 'pes', 'pmf', 'corr')))

    if mode in mode_list:
        #  GENERAL  ---------------------------------------------------
        h+=" -n | --name  <str>                basename for files\n"
        p.add_argument('-n', '--name', type=str)
        h+=" -c | --coord  <.crd>              coordinates file / path to files\n"
        p.add_argument('-c', '--coord', type=str)
        h+=" -s | --sys  <str>                 system basename to look for BIN/SELE\n"
        p.add_argument('-s', '--sys', type=str)
        h+=" --bin  <.bin>                     binary file for system\n"
        p.add_argument('--bin', type=str)
        h+=" --sele  <.dynn>                   selection file for QM and NOFIX\n\n"
        p.add_argument('--sele', type=str)
        
        if mode in ('scan', 'pes', 'irc', 'corr'):
            h+=" --out  <.out>                     file with distances and energies\n\n"
            p.add_argument('--out', type=str)

        #  COMPUTING  -------------------------------------------------
        h+=" --dry                             dry run: prepare but do not submit\n"
        p.add_argument('--dry', action='store_true')
        h+=" --queue  <str>                    queue to submit\n"
        p.add_argument('--queue', type=str)
        h+=" --cores  <int>                    total number of CPU cores\n"
        p.add_argument('--cores', type=int)
        h+=" --memory  <str>                   total amount of RAM memory\n\n"
        p.add_argument('--memory', type=str)

        #  QM  --------------------------------------------------------
        h+=" --charge  <int>                   QM-region charge\n"
        p.add_argument('--charge', type=int)
        h+=" --multi  <int>                    QM-region multiplicity\n"
        p.add_argument('--multi', type=int)
        h+=" --force_uhf                       force an unrestricted calculation\n\n"
        p.add_argument('--force_uhf', action='store_true')

        #  METHOD  ----------------------------------------------------
        h+=" --semiemp  <str>                  semi-empirical method\n"
        p.add_argument('--semiemp', type=str, choices=['AM1', 'RM1', 'PM3', 'MNDO', 'PDDG'])
        h+=" --gauss                           use Gaussian 09 software for the QM-region\n"
        p.add_argument('--gauss', action='store_true')
        h+=" --func  <str>                     DFT functional (for Gaussian)\n"
        p.add_argument('--func', type=str)
        h+=" --basis  <str>                    DFT basis set (for Gaussian)\n\n"
        p.add_argument('--basis', type=str)

    #  OPTIMIZATION  --------------------------------------------------
    if mode in ('mini', 'scan', 'pes'):
        h+=" {}:\n".format(mode.upper())
        h+=" --cg_steps  <int>                 Conjugate-Gradient maximum number of steps\n"
        p.add_argument('--cg_steps', type=int)
        h+=" --cg_tolerance  <float>           Conjugate-Gradient convergence criteria\n"
        p.add_argument('--cg_tolerance', type=float)
        h+=" --lbfgsb_steps  <int>             L-BFGSB maximum number of steps\n"
        p.add_argument('--lbfgsb_steps', type=int)
        h+=" --lbfgsb_tolerance  <float>       L-BFGSB convergence criteria\n\n"
        p.add_argument('--lbfgsb_tolerance', type=float)

    #  LOCATE  --------------------------------------------------------
    elif mode in ('locate'):
        h+=" {}:\n".format(mode.upper())
        h+=" --ts                              search a Transition State\n"
        p.add_argument('--ts', action='store_true')
        h+=" --loc_steps  <int>                baker's search maximum number of steps\n"
        p.add_argument('--loc_steps', type=int)
        h+=" --loc_tolerance  <float>          convergence criteria\n\n"
        p.add_argument('--loc_tolerance', type=float)

    #  IRC  -----------------------------------------------------------
    elif mode in ('irc'):
        h+=" {}:\n".format(mode.upper())
        h+=" --irc_dir  {-1,0,1}               initial direction to follow or both (0)\n"
        p.add_argument('--irc_dir', type=int, choices=[-1,0,1])
        h+=" --irc_steps  <int>                maximum number of steps\n"
        p.add_argument('--irc_steps', type=int)
        h+=" --irc_dsp  <float>                displacement on every step\n\n"
        p.add_argument('--irc_dsp', type=float)

    #  INTERACTION  ---------------------------------------------------
    elif mode in ('interaction'):
        h+=" {}:\n".format(mode.upper())
        h+=" --dcd_stride  <int>               read only every n-th frame of the trajectory\n"
        p.add_argument('--dcd_stride', type=int)
        h+=" --int_dcd  <.dcd>                 trajectory file along which calculate interactions\n\n"
        p.add_argument('--int_dcd', type=str)

    #  KIE  -----------------------------------------------------------
    elif mode in ('kie'):
        h+=" {}:\n".format(mode.upper())
        h+=" --temp  <float>                   temperature [k]\n"
        p.add_argument('--temp', type=float)
        h+=" --kie_atom  <str> <int> <str>     atom to calculate KIE (subsystem, residue number, atom name)\n"
        p.add_argument('--kie_atom', type=str)
        h+=" --kie_skip  <int>                 number of frequencies to skip\n"
        p.add_argument('--kie_skip', type=int)
        h+=" --kie_hess  <.dump>               hessian file\n\n"
        p.add_argument('--kie_hess', type=str)

    #  MD & PMF -------------------------------------------------------
    elif mode in ('md', 'pmf'):
        h+=" {}:\n".format(mode.upper())
        h+=" --temp  <float>                   temperature [k]\n"
        p.add_argument('--temp', type=float)
        h+=" --md_step  <float>                time step [ps]\n"
        p.add_argument('--md_step', type=float)
        h+=" --equi  <int>                     number of steps of equilibration\n"
        p.add_argument('--equi', type=int)
        h+=" --prod  <int>                     number of steps of production\n"
        p.add_argument('--prod', type=int)
        h+=" --dcd_freq  <int>                 frequency to save structures to the trajectory file\n"
        p.add_argument('--dcd_freq', type=int)
        h+=" --vel  <.vel>                     input velocities file to read instad of generate (for continuations)\n\n"
        p.add_argument('--vel', type=str)

    # dimensions
    if mode in ('pmf', 'corr'):
        h+=" -d | --dim  <int>                 dimension, otherwise guessed from config file\n\n"
        p.add_argument('-d', '--dim', type=int, default=0)

    # post-process
    if mode in ('scan', 'pes', 'pmf', 'corr'):
        h+=" -p | --post                       post-process routine after calculation\n\n"
        p.add_argument('-p', '--post', action='store_true')

    #  VERSION AND HELP  ----------------------------------------------
    p.add_argument('-v',
                   '--version',
                   action='version',
                   version='JACQUES v{} / GPL'.format(__version__))

    h+=" -v | --version                    show version info and exit\n"
    h+=" -h | --help                       display this help and exit\n"

    p.add_argument('-h',
                   '--help',
                   action='help',
                   default=argparse.SUPPRESS,
                   help=h)

    # return the built parser object
    return p
