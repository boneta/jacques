"""
=======================================================================
  File processing and post-processing routines
=======================================================================

  Functions
  ---------

    post_log
    post_out
    post_crd
    post_dat
    post

"""

import glob
import os
import re
import shutil
import sys

from .dynnconfig import DynnConfig
from .outfile import OutFile


def _natural_sort(l:list) -> list:
    '''Sort a list by natural order'''
    alphanum_key = lambda key: [int(c) if c.isdigit() else c.lower() for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)

def _progress_bar(name:str, n:int, n_tot:int, bar_length:int=21) -> None:
    '''Display a progress bar at a given step'''
    sys.stdout.write("\r# {}s  [{:21s}] - {:>6.2f}%".format(name,
                                                            "■"*int(bar_length*(n+1)/n_tot),
                                                            (n+1)/n_tot*100))
    sys.stdout.flush()

def post_log(files:list, file_final:str, rm:bool=True, progress_bar:bool=True) -> None:
    '''
        Post-Processing routine for .log files

        Parameters
        ----------
        files : list
            list of .log files
        file_final : str
            name of sorted .log file to write
        rm : bool, optional
            remove files after processing (def: True)
        progress_bar : bool, optional
            display progress bar (def: True)
    '''
    if file_final in files:
        files.remove(file_final)  # remove final file from file lists
    f_final = open(file_final, 'a+')
    for n, f in enumerate(files):
        f_final.write(f"\n#######  {f}  \n")
        with open(f, 'r') as f_tmp:
            f_final.write(f_tmp.read())
        if rm:
            os.remove(f)
        if progress_bar:
            _progress_bar("LOG", n, len(files))
    f_final.close()

def post_out(files:list, file_final:str=None, dynnconfig:'DynnConfig'=None, rm:bool=True) -> list:
    '''
        Post-Processing routine for .out files

        Parameters
        ----------
        files : list
            list of .out files
        file_final : str, optional
            name of sorted .out file to write
        dynnconfig : DynnConfig, optional
            DynnConfig object to read extreme indexes
        rm : bool, optional
            remove files after processing (def: True)
            only if file_final is not None and no missing values

        Returns
        -------
        ndarray of ndarray
            list of missing indexes
    '''
    out = OutFile(files)
    max_ndx = None
    if dynnconfig:
        dynnconfig.resolve_constr(out.dim)
        max_ndx = [dynnconfig.constr[i]['n'] for i in range(out.dim)]
    # remove duplicate index and find missing values
    out.unique()
    missing_ndx = out.missing_ndx(max_ndx)
    # write sorted file
    if file_final is not None:
        out.sort()
        out.zero()
        out.write(file_final)
    # remove files
    if rm and file_final is not None and len(missing_ndx) == 0:
        for f in files:
            os.remove(f)
    return missing_ndx

def post_crd(files:list, folder:str='crd', progress_bar:bool=True) -> None:
    '''
        Post-Processing routine for .crd files

        Parameters
        ----------
        files : list
            list of .crd files
        folder : str, optional
            name of folder to move .crd files to (def: 'crd')
        progress_bar : bool, optional
            display progress bar (def: True)
    '''
    if not os.path.isdir(folder):
        os.mkdir(folder)
    for n, f in enumerate(files):
        shutil.move(f, os.path.join(folder, f))
        if progress_bar:
            _progress_bar("CRD", n, len(files))

def post_dat(files:list, folder:str='dat', progress_bar:bool=True) -> None:
    '''
        Post-Processing routine for .dat files

        Parameters
        ----------
        files : list
            list of .dat files
        folder : str, optional
            name of folder to move .dat files to (def: 'dat')
        progress_bar : bool, optional
            display progress bar (def: True)
    '''
    if not os.path.isdir(folder):
        os.mkdir(folder)
    for n, f in enumerate(files):
        shutil.move(f, os.path.join(folder, f))
        if progress_bar:
            _progress_bar("DAT", n, len(files))

def post(dynnconfig:'DynnConfig', rm:bool=True) -> None:
    '''
        Post-Processing routine after a DYNAMON calculation

        Parameters
        ----------
        dynnconfig : DynnConfig
            DynnConfig object
        rm : bool, optional
            remove files after processing (def: True)
    '''

    mode_filetypes = {'scan': ('crd',),
                      'pes': ('log', 'out', 'crd'),
                      'pmf': ('log', 'out', 'crd', 'dat')}

    mode = dynnconfig.mode
    name = dynnconfig.name

    # check fundamental parameters
    if mode is None:
        sys.exit("ERROR: Missing MODE")
    elif mode not in mode_filetypes.keys():
        sys.exit(f"WARNING: No post-process routine for this mode '{mode}'. Nothing to do.")

    # list of files to process in natural order
    filetypes = mode_filetypes[mode]
    dir_files = dict()
    for filetype in filetypes:
        if filetype in ('log', 'crd', 'out'):
            file_pattern = f"{name}.*.{filetype}"
        elif filetype == 'dat':
            file_pattern = f"{filetype}_*"
        dir_files[filetype] = _natural_sort(glob.glob(file_pattern))

    sys.stdout.write(f"## POST-PROCESS: {mode}\n# NAME: {name}\n\n")

    for filetype in filetypes:
        # check no files matching
        if not dir_files[filetype]:
            sys.stdout.write(f"# {filetype.upper()}s  ⨯\n")
            continue
        # pre-loop actions
        if filetype == 'log':
            post_log(dir_files[filetype], f"{name}.log", rm)
        elif filetype == 'out':
            missing_out = post_out(dir_files[filetype], f"{name}.out", dynnconfig, rm=rm)
            if len(missing_out) > 0:
                sys.stdout.write(f"\r# OUTs ⨯  -  Missing {len(missing_out)} points\n")
                continue
        elif filetype == 'crd':
            post_crd(dir_files[filetype], folder='crd')
        elif filetype == 'dat':
            post_crd(dir_files[filetype], folder='dat')
        sys.stdout.write(f"\r# {filetype.upper()}s  ✔ {' ':<40}\n")

    sys.stdout.write("\n")
