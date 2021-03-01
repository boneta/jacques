#!/usr/bin/env python3

# Description: Convert dihedral parameters from AMBER to OPLS-AA
# Last update: 26-02-2021

import sys
import math as m
import argparse

##  DIHEDRAL PARAMETERS CLASS  ########################################
class DihedralParam:
    """
        Dihedral Angle - Parameters Class

        Attributes
        ----------
        e_thr : float
        angle_test : float
        phase_units : str
        p : dict
            dictionary of parameters
            'amber': {'dividier':[]*5,'barrier':[]*5,'phase':[]*5,'pn':[]*5}
            'opls' : []*5

        Methods
        -------
        __init__()
            Initialization of empty parameters
        list2amber(plist, phase_units)
            Read an input list of AMBER parameters
        list2opls(plist)
            Read an input list of OPLS parameters
        torsion_energy(ff, angle)
            Torsion energy value for an angle
        check_equal(angle_test, e_thr)
            Check if AMBER and OPLS parameters give equal torsion energies
        amber2opls()
            Convert AMBER parameters to OPLS
    """

    e_thr = 10e-6        # threshold value to consider energy equal
    angle_test = 2.1     # arbitrary angle to test fitting

    def __init__(self):
        """Initialization of empty parameters"""
        self.p = { 'amber' : { 'divider' : [0.0]*5,        # divider
                               'barrier' : [0.0]*5,        # barrier (named 'pk' / 'kd')
                               'phase'   : [0.0]*5,        # phase (rad)
                               'pn'      : [0]*5 } ,       # periodicity
                   'opls' : [ 0.0,                         # V0
                              0.0,                         # V1
                              0.0,                         # V2
                              0.0,                         # V3
                              0.0 ] }                      # V4

    def list2amber(self, plist, phase_units='degree') -> None:
        """
            Read an new input list of AMBER parameters

            Parameters
            ----------
            plist : list
                list of AMBER parameters, ordered and in sets of 4
                [IDIVF  PK  PHASE  PN ...]
            phase_units : {'rad', 'degree'}, optional
                angular units of the phase (def: degree)
        """
        self.__init__()
        self.phase_units = phase_units
        if len(plist) % 4 != 0: raise ValueError("Incomplete number of parameters (4 in every set)")
        for n in range( int(len(plist) / 4) ):
            pn = abs(int(plist[3+n*4]))         # periodicity number
            if not (0 <= pn <= 4):
                raise ValueError("Periodicity must be between 1 and 4")
            self.p['amber']['divider'][pn] = plist[0+n*4]
            self.p['amber']['barrier'][pn] = plist[1+n*4]
            self.p['amber']['phase'][pn] = plist[2+n*4]
            self.p['amber']['pn'][pn] = pn

    def list2opls(self, plist) -> None:
        """
            Read an input list of OPLS parameters

            Parameters
            ----------
            plist : list
                list of ordered OPLS parameters (can be missing from last)
                [V0  V1  V2  V3  V4]
        """
        self.__init__()
        if not (0 <= len(plist) < 5):
            raise ValueError("Number of OPLS parameters must be between 1 and 5")
        for i, param in enumerate(plist):
            self.p['opls'][i] = param

    def torsion_energy(self, ff, angle) -> float:
        """
            Torsion energy value for an angle

            Parameters
            ----------
            ff : {'gaff', 'amber', 'opls'}
                force field's function and parameters to use
            angle : float
                angle value to calculate energy [rad]

            Returns
            -------
            float
                energy value [kJ]
        """
        energy = 0
        # AMBER torsion energy
        if ff.lower() in ('gaff', 'amber'):
            param = self.p['amber']
            for i in range(len(param['divider'])) :
                if param['divider'][i] != 0 :
                    energy = energy + ( param['barrier'][i]/param['divider'][i]
                             * (1 +  m.cos(i*angle-param['phase'][i])) )
        # OPLS-AA torsion energy
        elif ff.lower() == 'opls':
            param = self.p['opls']
            energy = param[0]/2 + \
                     param[1]/2 * ( 1 + m.cos(1*angle) ) + \
                     param[2]/2 * ( 1 - m.cos(2*angle) ) + \
                     param[3]/2 * ( 1 + m.cos(3*angle) ) + \
                     param[4]/2 * ( 1 - m.cos(4*angle) )
        else:
            raise ValueError(f"Unkown force field: {ff}")
        return energy

    def check_equal(self, angle_test=angle_test, e_thr=e_thr) -> bool:
        """
            Check if AMBER and OPLS parameters give equal torsion energies

            Parameters
            ----------
            angle_test : float, optional
                arbitrary angle to test fitting [rad]
            e_thr : float, optional
                threshold value to consider energy equal

            Returns
            -------
            bool
                True if all tested angles are equal
        """
        angles = [0., m.pi/4, m.pi/2, m.pi, angle_test]
        e_diff = [abs(self.torsion_energy('amber',a)-self.torsion_energy('opls',a)) for a in angles]
        return all(e_diff) < e_thr

    def amber2opls(self) -> bool:
        """
            Convert AMBER parameters to OPLS

            Returns
            -------
            bool
                True if succeeded
        """
        opls  = self.p['opls']
        amber = self.p['amber']
        if self.phase_units == 'degree':
            amber['phase'] = [ m.radians(i) for i in amber['phase'] ]
        # Kasia style (Excel)
        for v in range(len(opls)) :
            if amber['divider'][v] != 0 :
                if v==1 and amber['phase'][v] == m.pi :
                    opls[v] = -2* amber['barrier'][v] / amber['divider'][v]
                else:
                    opls[v] = 2* amber['barrier'][v] / amber['divider'][v]
        opls[0] = 2* ( - self.torsion_energy('opls', 0.) + self.torsion_energy('amber', 0.) )
        # Alternative method (invert and correct 0)
        if not self.check_equal():
            opls[0] = 0
            opls = [ -1*i if i !=0 else 0 for i in opls ]
            opls[0] = 2* ( - self.torsion_energy('opls', 0.) + self.torsion_energy('amber', 0.) )
        # return success status
        return self.check_equal()


##  MAIN  #############################################################
if __name__ == '__main__':

    # parser
    parser = argparse.ArgumentParser(description='Convert AMBER parameters to OPLS format\n\n'+
                                                 '  AMBER Input:    IDIVF  PK  PHASE  PN\n'+
                                                 '  OPLS Output:    V0  V1  V2  V3  V4',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('amber_param',
                        metavar='X X X X', type=float, nargs='+',
                        help='Amber parameters (IDIVF PK PHASE PN)')
    parser.add_argument('--phase', metavar='rad | degree',
                        choices=['rad', 'degree'], default='degree',
                        help='rad or degree input (def: degree)')
    args = parser.parse_args()

    # check number of parameters
    if len(args.amber_param) % 4 != 0:
        sys.exit("ERROR: Incomplete number of parameters (4 in every set)")

    # read and convert parameters
    dihedral = DihedralParam()
    dihedral.list2amber(args.amber_param, args.phase)
    dihedral.amber2opls()

    # print conversion
    if not dihedral.check_equal():
        sys.stdout.write("WARNING: Parameters do not match. \n WARNING: Last Parameters Calculated: ")
    sys.stdout.write("{:>8.3f} {:>8.3f} {:>8.3f} {:>8.3f} {:>8.3f} \n".format(*dihedral.p['opls']))
