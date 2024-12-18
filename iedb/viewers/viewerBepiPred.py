# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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
import os, webbrowser

from pyworkflow.protocol import params, Protocol

from pwchem.viewers.viewers_sequences import SequenceGeneralViewer

from ..protocols import ProtBepiPredPrediction

class ViewerBepiPred(SequenceGeneralViewer):
  _label = 'BepiPred viewer'
  _targets = [ProtBepiPredPrediction]

  def _defineParams(self, form):
    super()._defineParams(form)
    group = form.addGroup('BepiPred scores view')
    group.addParam('displayBepiPred', params.LabelParam,
                   label='Display with BepiPred scores: ',
                   help='Display html with raw BepiPred scores')

  def _getVisualizeDict(self):
    vDic = super()._getVisualizeDict()
    vDic.update({'displayBepiPred': self._showBepiPred})
    return vDic

  def _showBepiPred(self, paramName=None):
      htmlFile = self.protocol._getExtraPath('output_interactive_figures.html')
      webbrowser.open(htmlFile)

  def getOutDir(self):
    return os.path.abspath(self.getProtocol()._getExtraPath()
                           if self.getProtocol() else self.getProject().getTmpPath())

  def getProtocol(self):
    if hasattr(self, 'protocol') and isinstance(self.protocol, Protocol):
      return self.protocol

  def getOutputSet(self):
    return self.protocol.outputROIs