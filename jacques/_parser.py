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


def parser(mode=None) -> 'argparse.Namespace':
    '''
        Construct specific parser for a MODE

        Parameters
        ----------
        mode : str, optional
            calculation mode

        Returns
        -------
        argparse.Namespace
            parsed arguments
    '''

    mode_list = ('sp', 'mini', 'locate', 'md', 'interaction',
                 'kie', 'irc', 'pes', 'scan', 'pmf', 'corr')

    # custom argparse to include help
    class ArgumentParserCustom(argparse.ArgumentParser):
        def format_help(self):
             return h
    p = ArgumentParserCustom(formatter_class=argparse.RawTextHelpFormatter, add_help=False)

    # processed 'mode'
    mode = mode.lower() if mode.lower() in mode_list else '<mode>'

    # help text
    h=  " ##################  JACQUES  ##################\n"
    h+=" -- Friendly interface for QM/MM calculations --\n\n"
    h+=f" USAGE:   jacques {mode} [[--option arg] ...]\n\n"
    h+=" OPTIONS:\n"

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
        h+=" -c | --coord  <.crd>              coordinates file / path to files\n"
        h+=" -s | --sys  <str>                 system basename to look for BIN/SELE\n"
        h+=" --bin  <.bin>                     binary file for system\n"
        h+=" --sele  <.dynn>                   selection file for QM and NOFIX\n\n"
        if mode in ('scan', 'pes', 'irc', 'corr'):
            h+=" --out  <.out>                     file with distances and energies\n\n"
        #  COMPUTING  -------------------------------------------------
        h+=" --dry                             dry run: prepare but do not submit\n"
        h+=" --queue  <str>                    queue to submit\n"
        h+=" --cores  <int>                    total number of CPU cores\n"
        h+=" --memory  <str>                   total amount of RAM memory\n\n"
        #  QM  --------------------------------------------------------
        h+=" --charge  <int>                   QM-region charge\n"
        h+=" --multi  <int>                    QM-region multiplicity\n"
        h+=" --force_uhf                       force an unrestricted calculation\n\n"
        #  METHOD  ----------------------------------------------------
        h+=" --semiemp  <str>                  semi-empirical method\n"
        h+=" --gauss                           use Gaussian 09 software for the QM-region\n"
        h+=" --func  <str>                     DFT functional (for Gaussian)\n"
        h+=" --basis  <str>                    DFT basis set (for Gaussian)\n\n"
    p.add_argument('-n', '--name', type=str)
    p.add_argument('-c', '--coord', type=str)
    p.add_argument('-s', '--sys', type=str)
    p.add_argument('--bin', type=str)
    p.add_argument('--sele', type=str)
    p.add_argument('--out', type=str)
    p.add_argument('--dry', action='store_true')
    p.add_argument('--queue', type=str)
    p.add_argument('--cores', type=int)
    p.add_argument('--memory', type=str)
    p.add_argument('--charge', type=int)
    p.add_argument('--multi', type=int)
    p.add_argument('--force_uhf', action='store_true')
    p.add_argument('--semiemp', type=str, choices=['AM1', 'RM1', 'PM3', 'MNDO', 'PDDG'])
    p.add_argument('--gauss', action='store_true')
    p.add_argument('--func', type=str)
    p.add_argument('--basis', type=str)

    #  OPTIMIZATION  --------------------------------------------------
    if mode in ('mini', 'scan', 'pes'):
        h+=" {}:\n".format(mode.upper())
        h+=" --cg_steps  <int>                 Conjugate-Gradient maximum number of steps\n"
        h+=" --cg_tolerance  <float>           Conjugate-Gradient convergence criteria\n"
        h+=" --lbfgsb_steps  <int>             L-BFGSB maximum number of steps\n"
        h+=" --lbfgsb_tolerance  <float>       L-BFGSB convergence criteria\n\n"
    p.add_argument('--cg_steps', type=int)
    p.add_argument('--cg_tolerance', type=float)
    p.add_argument('--lbfgsb_steps', type=int)
    p.add_argument('--lbfgsb_tolerance', type=float)

    #  LOCATE  --------------------------------------------------------
    if mode in ('locate'):
        h+=" {}:\n".format(mode.upper())
        h+=" --ts                              search a Transition State\n"
        h+=" --loc_steps  <int>                baker's search maximum number of steps\n"
        h+=" --loc_tolerance  <float>          convergence criteria\n\n"
    p.add_argument('--ts', action='store_true')
    p.add_argument('--loc_steps', type=int)
    p.add_argument('--loc_tolerance', type=float)

    #  IRC  -----------------------------------------------------------
    if mode in ('irc'):
        h+=" {}:\n".format(mode.upper())
        h+=" --irc_dir  {-1,0,1}               initial direction to follow or both (0)\n"
        h+=" --irc_steps  <int>                maximum number of steps\n"
        h+=" --irc_dsp  <float>                displacement on every step\n\n"
        h+=" IRC POST-PROCESSING:\n"
        h+=" --irc_invert                      invert the direction of the IRC found\n"
        h+=" --irc_dat  <.dat>                 unified reaction profile file to write\n"
        h+=" --irc_crd  <str>                  extract the IRC coordinates as .crd in this folder ('coord' needed)\n"
        h+="                                   number of coordinates to extract is taken from 'irc_steps', use '-1' to extract all \n\n"
    p.add_argument('--irc_dir', type=int, choices=[-1,0,1])
    p.add_argument('--irc_steps', type=int)
    p.add_argument('--irc_dsp', type=float)
    p.add_argument('--irc_invert', action='store_true')
    p.add_argument('--irc_dat', type=str)
    p.add_argument('--irc_crd', type=str)

    #  INTERACTION  ---------------------------------------------------
    if mode in ('interaction'):
        h+=" {}:\n".format(mode.upper())
        h+=" --dcd_stride  <int>               read only every n-th frame of the trajectory\n"
        h+=" --int_dcd  <.dcd>                 trajectory file along which calculate interactions\n\n"
    p.add_argument('--dcd_stride', type=int)
    p.add_argument('--int_dcd', type=str)

    #  KIE  -----------------------------------------------------------
    if mode in ('kie'):
        h+=" {}:\n".format(mode.upper())
        h+=" --temp  <float>                   temperature [k]\n"
        h+=" --kie_atom  <str> <int> <str>     atom to calculate KIE (subsystem, residue number, atom name)\n"
        h+=" --kie_skip  <int>                 number of frequencies to skip\n"
        h+=" --kie_hess  <.dump>               hessian file\n\n"
    p.add_argument('--temp', type=float)
    p.add_argument('--kie_atom', type=str)
    p.add_argument('--kie_skip', type=int)
    p.add_argument('--kie_hess', type=str)

    #  MD & PMF -------------------------------------------------------
    if mode in ('md', 'pmf'):
        h+=" {}:\n".format(mode.upper())
        h+=" --temp  <float>                   temperature [k]\n"
        h+=" --md_step  <float>                time step [ps]\n"
        h+=" --equi  <int>                     number of steps of equilibration\n"
        h+=" --prod  <int>                     number of steps of production\n"
        h+=" --dcd_freq  <int>                 frequency to save structures to the trajectory file\n"
        h+=" --vel  <.vel>                     input velocities file to read instad of generate (for continuations)\n\n"
    p.add_argument('--md_step', type=float)
    p.add_argument('--equi', type=int)
    p.add_argument('--prod', type=int)
    p.add_argument('--dcd_freq', type=int)
    p.add_argument('--vel', type=str)

    # dimensions
    if mode in ('pmf', 'corr'):
        h+=" -d | --dim  <int>                 dimension, otherwise guessed from config file\n\n"
    p.add_argument('-d', '--dim', type=int, default=0)

    # post-process
    if mode in ('scan', 'pes', 'pmf', 'corr', 'irc'):
        h+=" -p | --post                       perform post-process routine for after calculation instead of launching\n\n"
    p.add_argument('-p', '--post', action='store_true')

    #  VERSION AND HELP  ----------------------------------------------
    h+=" -v | --version                    show version info and exit\n"
    h+=" -h | --help                       display this help and exit\n"
    p.add_argument('-v',
                   '--version',
                   action='version',
                   version='JACQUES v{} / GPL'.format(__version__))
    p.add_argument('-h',
                   '--help',
                   action='help',
                   default=argparse.SUPPRESS,
                   help=h)

    # return the arguments
    return p.parse_args()

