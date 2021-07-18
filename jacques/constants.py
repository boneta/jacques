"""
=======================================================================
  Constants & Conversions
=======================================================================

Thermodynamical constants and conversions
http://physics.nist.gov/cuu/Constants/index.html

Temperature conversion
https://www.nist.gov/pml/weights-and-measures/si-units-temperature


Functions
---------

    temp2temp

"""


##  Thermodynamical constants
C         = 299792458.                # m * s-1                 # speed of light in vacuum
H         = 6.62607015e-34            # J * s                   # Planck constant
KB        = 1.380649e-23              # J * K-1                 # Boltzmann constant
NA        = 6.02214076e23             # mol-1                   # Avogadro constant
DA        = 1.6605390660e-27          # kg                      # atomic mass constant
ME        = 9.1093837015e-31          # kg                      # electron mass
R         = KB * NA                   # J * K-1 * mol-1         # molar gas constant
BOHR      = 5.291772109e-11           # m                       # Bohr radius (a0)


##  Conversions
CAL2J     = 4.184                     # cal -> J
J2CAL     = 1./CAL2J                  # J -> cal
EV2J      = 1.602176634e-19           # eV -> J
HA2J      = 4.3597447222071e-18       # Hartree -> J
HA2JMOL   = HA2J * NA                 # Hartree -> J * mol-1
A2NM      = 0.1                       # Å -> nm
BOHR2A    = BOHR * 1e10               # a0 -> Å
A2BOHR    = 1./BOHR2A                 # Å -> a0


def temp2temp(temp, in_u, out_u='K'):
    """
        Convert temperature values among units {K, C, F}

        Parameters
        ----------
        temp : float
            temperature value to convert, in 'in_u' units
        in_u : {K, C, F}
            input units
        out_u : {K, C, F}, optional
            output units (def: K)

        Returns
        -------
        t_conv : float
            temperature value converted, in 'out_u' units
    """

    in_u  = in_u.upper()
    out_u = out_u.upper()

    # standarize to K
    if in_u == 'K':
        t_std = temp
    elif in_u == 'C':
        t_std = temp + 273.15
    elif in_u == 'F':
        t_std = ( temp - 32 ) / 1.8 + 273.15
    else:
        raise ValueError(f'Unrecognized temperature unit: {in_u}')

    # convert
    if out_u == 'K':
        t_conv = t_std
    elif out_u == 'C':
        t_conv = t_std - 273.15
    elif out_u == 'F':
        t_conv = ( t_std - 273.15 ) * 1.8 + 32
    else:
        raise ValueError(f'Unrecognized temperature unit: {out_u}')

    return t_conv
