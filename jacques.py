#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#
#
#   _     _     _            _
#  | | _ | | _ | |          (_)_____   _____ ____   __  __ ___ _______
#  )=(| | _ | |)=(         / /_     ` / ___// __ ` / / / /  _ \_  ___/
# (   )=(| |)=(   )       / / / <3  // /__ / /_/ // /_/ //  __/(__  )
#  \_(   )=(   )_/       / /  \__,_/ \___/ \__  / \____/ \___//____/
#     \_(   )_/      ___/ /                  /_/
#        \_/        /____/   Friendly interface for QM/MM calculations
#
#
#
#
#######################################################################
##                                                                   ##
##                              JACQUES                              ##
##                                                                   ##
#######################################################################
#
# Copyright (C) 2021, Sergio Boneta
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see https://www.gnu.org/licenses/
#

import os
import sys

from jacques          import jacques_path
from jacques.__main__ import main

if __name__ == '__main__':
    # check for tools call
    if len(sys.argv)>1 and sys.argv[1].lower() in ("tools", "t"):
        tools_path = jacques_path+'/../tools/'
        if len(sys.argv) == 2:
            sys.stdout.write("ERROR: No tool specified. Use -h for help.\n")
        elif sys.argv[2].lower() in ("-h", "--help"):
            with open(tools_path+'README.md', 'r') as helpfile:
                sys.stdout.writelines([line for line in helpfile.readlines() if "``" not in line])
        else:
            os.system(tools_path+' '.join(sys.argv[2:]))
    else:
        main()
