<img width="350" height="150" src="./docs/jacques_logo.svg" align="right" />

# JACQUES

![python](https://img.shields.io/badge/python-3-red.svg)
![Platform](https://img.shields.io/badge/platform-linux-lightgrey.svg)
[![License: GPLv3](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) \
*Friendly interface for QM/MM calculations*

## Installation
This installation also includes the [DYNAMO<sup>N</sup>](https://github.com/boneta/dynamon) package.

1. Clone the repository in your preferred location.
2. Set the environment variables of JACQUES (and DYNAMO<sup>N</sup>). For a better usability, *source* the file `jacques.rc` in your *.bashrc*.
3. Compile DYNAMON with a Fortran compiler. Tested with *gfortran*.

```bash
git clone --recursive https://github.com/boneta/jacques
source jacques/jacques.rc
make -C $DYNAMON/src
```

## Usage
For more information, see the [documentation](./docs/README.md).
```
jacques <mode> [options]
```
