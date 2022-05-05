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
import shutil
import sys

import numpy as np

from .dynnconfig import DynnConfig, _natural_sort


def _progress_bar(name, n, n_tot, bar_length=21) -> None:
    sys.stdout.write("\r# {}s  [{:21s}] - {:>6.2f}%".format(name,
                                                            "■"*int(bar_length*(n+1)/n_tot),
                                                            (n+1)/n_tot*100))
    sys.stdout.flush()

def post_log(files:list, file_final, rm=True, progress_bar=True) -> None:
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

def post_out(files, file_final=None, dynnconfig=None, ndx=None, rm=True) -> list:
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
        ndx : list, optional
            list of extreme indexes to find
            default behavior is from 0 to max index found
        rm : bool, optional
            remove files after processing (def: True)
            only if file_final is not None and no missing values

        Returns
        -------
        list
            list of missing points
    '''
    # read all files
    out_data_all = []
    for out in files:
        with open(out, 'r') as f:
            out_first_line = f.readline().strip()
            out_data_all.append(np.loadtxt(f, dtype=float, comments='#'))
    out_data = np.concatenate(out_data_all)

    # guess dimension (number of columns)
    if out_data.shape[1] == 5:
        dim = 1
    elif out_data.shape[1] == 8:
        dim = 2
    else:
        raise ValueError("Inconsistent number of columns in .out file")

    if dynnconfig:
        dynnconfig.resolve_constr(dim)

    # remove duplicate index, sort data and find missing values
    if dim == 1:
        out_data = out_data[np.unique(out_data[:, 3], return_index=True)[1]]
        out_data = out_data[out_data[:, 3].argsort()]
        if ndx is not None:
            ndx = (0, ndx[0])
        elif dynnconfig:
            ndx = (0, dynnconfig.constr[0]['n'])
        else:
            ndx = (0, max(out_data[:, 3]))
        total_values = set(np.arange(ndx[0], ndx[1]+1))
        missing_values = total_values - set(out_data[:, 3])
        missing_values = [str(int(i)) for i in missing_values]
        missing_values.sort()
    elif dim == 2:
        out_data = out_data[np.unique(out_data[:, 4:6], axis=0, return_index=True)[1]]
        out_data = out_data[np.lexsort((out_data[:, 5], out_data[:, 4]))]
        if ndx is not None and len(ndx) >= 2:
            ndx = ((0, ndx[0]), (0, ndx[1]))
        elif dynnconfig:
            ndx = ((0, dynnconfig.constr[0]['n']), (0, dynnconfig.constr[1]['n']))
        else:
            ndx = ((0, max(out_data[:, 4])), (0, max(out_data[:, 5])))
        total_values = {(i,j) for i in np.arange(ndx[0][0], ndx[0][1]+1) for j in np.arange(ndx[1][0], ndx[1][1]+1)}
        missing_values = total_values - {(i,j) for i,j in zip(out_data[:, 4], out_data[:, 5])}
        missing_values = [str(int(i))+' '+str(int(j)) for i,j in missing_values]
        missing_values.sort()

    if file_final is not None:
        # minimum energies to zero
        for i in [dim+0, dim+1]:
            out_data[:, i] -= min(out_data[:, i])
        # save sorted output file
        fmt = r'%12.4f  '*dim + r'%20.10f  %20.10f  ' + r'%5d  '*dim + r'%12.4f  '*dim
        np.savetxt(file_final, out_data, fmt=fmt, comments='#', header=out_first_line[1:])

    if rm and file_final is not None and len(missing_values) == 0:
        for out in files:
            os.remove(out)

    return missing_values

def post_crd(files:list, folder='crd', progress_bar=True) -> None:
    if not os.path.isdir(folder):
        os.mkdir(folder)
    for n, f in enumerate(files):
        shutil.move(f, os.path.join(folder, f))
        if progress_bar:
            _progress_bar("CRD", n, len(files))

def post_dat(files:list, folder='dat', progress_bar=True) -> None:
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

    mode_filetypes = {'scan': ('crd'),
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
            if missing_out:
                sys.stdout.write(f"\r# OUTs ⨯  -  Missing {len(missing_out)} points\n")
                continue
        elif filetype == 'crd':
            post_crd(dir_files[filetype], folder='crd')
        elif filetype == 'dat':
            post_crd(dir_files[filetype], folder='dat')
        sys.stdout.write(f"\r# {filetype.upper()}s  ✔ {' ':<40}\n")

    sys.stdout.write("\n")
