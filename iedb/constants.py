# coding: latin-1
# **************************************************************************
# *
# * Authors:  Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *
# * Biocomputing Unit, CNB-CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

# Common constants
DEFAULT_VERSION = '3.0'

# Package dictionaries
BEPIPRED_DIC =  {'name': 'BepiPred',    'version': '3.0', 'pattern': 'bepipred',
                 'home': 'BEPIPRED_HOME', 'zip': 'BEPIPRED_ZIP', 'activation': 'BEPIPRED_ACTIVATION_CMD'}

MHCI_DIC = {'name': 'mhc_i',    'version': '3.1.5', 'pattern': 'mhc_i',
            'home': 'MHC_I_HOME', 'tar': 'MHC_I_TAR'}

MHCII_DIC = {'name': 'mhc_ii',    'version': '3.1.11', 'pattern': 'mhc_ii',
             'home': 'MHC_II_HOME', 'tar': 'MHC_II_TAR'}

COVE_DIC = {'name': 'population_coverage',    'version': '3.0.2', 'pattern': 'population_coverage',
            'home': 'COVERAGE_HOME', 'tar': 'COVERAGE_TAR'}

READ_URL = 'https://github.com/scipion-chem/scipion-chem-bepipred'

NOINSTALL_WARNING = f'Installation could not be completed because BepiPred download or installation has not been found.\n' \
                    f'Please check the scipion-chem-bepipred README file to see more details about how to proceed with ' \
                    f'the installation. You can find this instruction in {READ_URL}'
