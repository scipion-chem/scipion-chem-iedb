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

from pwchem.viewers import SequenceGeneralViewer, ViewerGeneralStructROIs
from pwchem.objects import SetOfSequenceROIs

from ..protocols import ProtElliProPrediction

class ViewerElliPro(SequenceGeneralViewer, ViewerGeneralStructROIs):
  _label = 'ElliPro viewer'
  _targets = [ProtElliProPrediction]
  _seqOutput = [SetOfSequenceROIs]

  def __init__(self, **kwargs):
    super().__init__(**kwargs)

  def _getVisualizeDict(self):
    vDic = SequenceGeneralViewer._getVisualizeDict(self)
    vDic.update(ViewerGeneralStructROIs._getVisualizeDict(self))
    return vDic

  def _defineParams(self, form):
    SequenceGeneralViewer._defineParams(self, form)
    ViewerGeneralStructROIs._defineParams(self, form)

  def getOutSequences(self):
    if self.checkIfProtocol():
      for oAttr in self.protocol.iterOutputAttributes():
        for oType in self._seqOutput:
          if isinstance(getattr(self.protocol, oAttr[0]), oType):
            return getattr(self.protocol, oAttr[0])
    else:
      return self.protocol

