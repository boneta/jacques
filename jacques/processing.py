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
    post_irc
    post

"""

import glob
import os
import re
import shutil
import sys

import numpy as np
from pdb4all import PDB, Traj

from .dynnconfig import DynnConfig
from .outfile import OutFile


def _natural_sort(l:list) -> list:
    '''Sort a list by natural order'''
    alphanum_key = lambda key: [int(c) if c.isdigit() else c.lower() for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)

def _progress_bar(name:str, n:int, n_tot:int, bar_length:int=21) -> None:
    '''Display a progress bar at a given step'''
    sys.stdout.write("\r# {}   [{:21s}] - {:>6.2f}%".format(name,
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

def post_irc(name:str, invert:bool=False, irc_dat:str=None, irc_crd:str=None, coord:str=None, num_coord:int=100, progress_bar:bool=True) -> None:
    '''
        Post-Processing routine for IRC calculations

        Parameters
        ----------
        name : str
            name of IRC calculation (used to look up folders and files)
        invert : bool, optional
            invert the direction of the IRC found (def: False)
        irc_dat : str, optional
            name of unified reaction profile file to write (.dat)
        irc_crd : str, optional
            folder to extract IRC coordinates to as .crd files
        coord : str, optional
            coordinates file to get topology from (.crd)
        num_coord : int, optional
            number of coordinates to extract on each side, geometrically spaced from the TS (def: 100)
            use -1 to extract all coordinates
        progress_bar : bool, optional
            display progress bar (def: True)
    '''
    def geomrange(start:int, stop:int, num:int) -> np.ndarray:
        '''Return n unique integer values geometrically spaced within a given interval'''
        if num > stop-start+1:
            raise ValueError("num must be smaller than stop-start+1")
        if start == 0:
            start = 1e-6
        l = np.array([])
        num_geom = num
        while len(l) < num:
            l = np.unique(np.geomspace(start, stop, num=num_geom, dtype=int))
            num_geom += 1
        return l

    # define folders/files
    irc_for = f"{name}-FOR"
    irc_back = f"{name}-BACK"
    for f in (irc_for, irc_back):
        if not os.path.isdir(f):
            sys.exit(f"ERROR: Folder '{f}' not found")

    # unify reaction profile (.dat)
    if irc_dat:
        sys.stdout.write(f"\r# Unified reaction profile:  {irc_dat}\n")
        dat_files = [f"{irc_back}/{irc_back}.dat", f"{irc_for}/{irc_for}.dat"]
        dat_data = []
        for f in dat_files:
            if not os.path.isfile(f):
                sys.exit(f"ERROR: File '{f}' not found")
            dat_data.append(np.loadtxt(f))
        dat_data = np.dstack(dat_data).swapaxes(0, 2)
        dat_data[:,0,:] += 1        # start index at 1
        if invert:
            dat_data = np.flip(dat_data, axis=0)
        dat_data[0] = np.flip(dat_data[0], axis=1)  # reverse the back
        dat_data[0,0,:] *= -1                       # negative index for back
        dat_data = np.hstack(dat_data).T
        np.savetxt(irc_dat, dat_data, fmt=' %6d    %f')

    # extract IRC coordinates (.crd)
    if irc_crd and coord:
        sys.stdout.write(f"\r# Extracting IRC coordinates to '{irc_crd}'\n")
        dcd_files = [f"{irc_back}/{irc_back}.dcd", f"{irc_for}/{irc_for}.dcd"]
        if invert:
            dcd_files = dcd_files[::-1]
        os.makedirs(irc_crd, exist_ok=True)
        ncrd = 1
        traj = Traj()
        for nside, dcd_file in enumerate(dcd_files, 1):
            traj.read_dcd(dcd_file, coord)
            nframes = traj.nframes
            num_coord_side = nframes if num_coord == -1 else min(num_coord, nframes)
            traj.frames_xyz = [traj.frames_xyz[i] for i in geomrange(0, nframes-1, num_coord_side)]
            if nside == 1:
                traj.frames_xyz = traj.frames_xyz[::-1]
            for i, frame in zip(range(ncrd, num_coord_side+ncrd), traj):
                frame.write_crd(f"{irc_crd}/{name}-{i:d}.crd")
                if progress_bar:
                    _progress_bar(f"CRD-{nside}", i-ncrd, num_coord_side)
            ncrd += num_coord_side
        sys.stdout.write(f"\r{' ':<80}\n\n")

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
                      'pmf': ('log', 'out', 'crd', 'dat'),
                      'irc': None}

    mode = dynnconfig.mode
    name = dynnconfig.name

    # check fundamental parameters
    if mode is None:
        sys.exit("ERROR: Missing MODE")
    elif mode not in mode_filetypes.keys():
        sys.exit(f"WARNING: No post-process routine for this mode '{mode}'. Nothing to do.")

    sys.stdout.write(f"## POST-PROCESS: {mode}\n# NAME: {name}\n\n")

    # irc post-processing
    if mode == 'irc':
        post_irc(name,
                 dynnconfig.opt['irc_invert'],
                 dynnconfig.opt['irc_dat'],
                 dynnconfig.opt['irc_crd'],
                 dynnconfig.opt['coord'],
                 dynnconfig.opt['irc_ncrd'])
        return

    # list of files to process in natural order
    filetypes = mode_filetypes[mode]
    dir_files = dict()
    for filetype in filetypes:
        if filetype in ('log', 'crd', 'out'):
            file_pattern = f"{name}.*.{filetype}"
        elif filetype == 'dat':
            file_pattern = f"{filetype}_*"
        dir_files[filetype] = _natural_sort(glob.glob(file_pattern))

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
