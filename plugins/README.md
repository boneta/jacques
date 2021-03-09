# Plug-ins

## [PyMOL](https://pymolwiki.org) `dynamon_pymol.py`

The PyMOL plug-in adds capabilities to this molecular viewer to manage files in fDynamo/DYNAMO format.

#### Installation

It can be used as a generic script/plug-in, but for an easiest usage it's advised to install it through the [plugin manager](https://pymolwiki.org/index.php/Plugin_Manager):

`Plugin > Plugin Manager > Install New Plugin > Install from local file`

#### Usage

  - Write QM atoms and NOFIX residues from selection: \
        `write_qm  dynn_file [, selection ]` \
        `write_nofix  dynn_file [, selection ]`
