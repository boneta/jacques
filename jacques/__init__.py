
r'''.
   _     _     _            _
  | | _ | | _ | |          (_)_____   _____ ____   __  __ ___ _______
  )=(| | _ | |)=(         / /_     ` / ___// __ ` / / / /  _ \_  ___/
 (   )=(| |)=(   )       / / / <3  // /__ / /_/ // /_/ //  __/(__  )
  \_(   )=(   )_/       / /  \__,_/ \___/ \__  / \____/ \___//____/
     \_(   )_/      ___/ /                  /_/
        \_/        /____/   Friendly interface for QM/MM calculations

'''

import os
import sys
import json


__version__ = '0.1.0'


##  initial configuration  ############################################

# absolute path of JACQUES package
jacques_path = os.path.dirname(os.path.realpath(__file__))

# absolute path of DYNAMON
dynamon_path = os.path.expandvars('$DYNAMON')

# read settings
try:
    with open(jacques_path+'/../settings.json', 'r') as _file:
        settings = json.loads(_file.read())
except:
    sys.stdout.write("ERROR: Problems reading 'settings.json' file \n")
    sys.exit(1)

# default queues names
queues_def = settings['queues']

# default options
opt_def = settings['options']

# default paths to store jobs and queue messages
msgf = os.path.expandvars(settings['folders']['messages'])
jobf = os.path.expandvars(settings['folders']['jobs'])
for _folder in msgf, jobf:
    if not os.path.exists(_folder): os.makedirs(_folder)
