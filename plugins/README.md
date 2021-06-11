# Plug-ins

## **[PyMOL](https://pymolwiki.org)**

`dynamon_pymol.py` The PyMOL plug-in adds capabilities to this molecular viewer to manage files in fDynamo/DYNAMO format.

### Installation

It can be used as a generic script/plug-in, but for an easiest usage it's advised to install it through the [plugin manager](https://pymolwiki.org/index.php/Plugin_Manager): `Plugin > Plugin Manager > Install New Plugin`

### Usage

  - Write a selection of atoms to a file: \
        `write_sele filename [, selection_name [, selection [, resolution ]]]`

  - Write directly QM atoms or NOFIX residues from selection: \
        `write_qm  filename [, selection ]` \
        `write_nofix  filename [, selection ]`

  - Load/read extended file formats:  `load  filename [...]`
    - .crd   -  fDynamo's coordinates file
    - .dynn  -  DYNAMON options/selection file
