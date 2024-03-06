# -*- coding: utf-8 -*-
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

"""
"""

from pyworkflow.gui import ListTreeProviderString, dialog
from pyworkflow.object import String

from pwchem.wizards import VariableWizard

from ..protocols import ProtMHCIPrediction

#################### MHC wizards #####################

class SelectMultiAlleleWizard(VariableWizard):
  _targets, _inputs, _outputs = [], {}, {}

  def getListOfElements(self, form):
    return form.protocol.getAvailableAlleles()

  def show(self, form, *params):
    inputParam, outputParam = self.getInputOutput(form)
    elementList = self.getListOfElements(form)

    finalList = []
    for i in elementList:
      finalList.append(String(i))
    provider = ListTreeProviderString(finalList)
    dlg = dialog.ListDialog(form.root, "Allele for prediction", provider,
                            "Select one or several alleles")

    vals = [v.get() for v in dlg.values]
    form.setVar(outputParam[0], ', '.join(vals))

SelectMultiAlleleWizard().addTarget(protocol=ProtMHCIPrediction,
                                    targets=['alleleCustom'],
                                    inputs=[],
                                    outputs=['alleleCustom'])

class SelectMultiLengthWizard(VariableWizard):
  _targets, _inputs, _outputs = [], {}, {}

  def getListOfElements(self, form):
    return list(map(str, range(8, 15)))

  def show(self, form, *params):
    inputParam, outputParam = self.getInputOutput(form)
    elementList = self.getListOfElements(form)

    finalList = []
    for i in elementList:
      finalList.append(String(i))
    provider = ListTreeProviderString(finalList)
    dlg = dialog.ListDialog(form.root, "Epitope lengths for prediction", provider,
                            "Select one or several lengths")

    vals = [v.get() for v in dlg.values]
    form.setVar(outputParam[0], ', '.join(vals))

SelectMultiLengthWizard().addTarget(protocol=ProtMHCIPrediction,
                                    targets=['lengths'],
                                    inputs=[],
                                    outputs=['lengths'])
