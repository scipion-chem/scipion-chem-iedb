# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

import os
import numpy as np

from pwem.protocols import EMProtocol
from pyworkflow.protocol import params

from pwchem.objects import Sequence, SequenceROI, SetOfSequenceROIs

from .. import Plugin as iedbPlugin
from ..constants import POP_DIC
from ..utils import getAllMHCIIAlleles

class ProtMHCIIPopulationCoverage(EMProtocol):
  """Calculates the population coverage of a series of MHC epitopes"""
  _label = 'mhc population coverage'


  def __init__(self, **kwargs):
    EMProtocol.__init__(self, **kwargs)

  def _defineParams(self, form):
    form.addSection(label='Input')
    iGroup = form.addGroup('Input')
    iGroup.addParam('inputSequenceROIs', params.PointerParam, pointerClass="SetOfSequenceROIs",
                    label='Input MHC epitopes: ',
                    help="Set of sequence MHC epitopes to calculate the population coverage on")

    iGroup.addParam('mhc', params.EnumParam, label='MHC type: ', display=params.EnumParam.DISPLAY_HLIST,
                    default=2, choices=['I', 'II', 'combined'],
                    help="Host specie  to predict the MHC-II epitopes on.")

    pGroup = form.addGroup('Population')
    pGroup.addParam('pop', params.EnumParam, label='Population by: ',
                    default=0, choices=['Area', 'Ethnicity'],
                    help="Choose population by geographic area or ethnicity")

    pGroup.addParam('area1', params.StringParam, label='Continent: ', default='All', condition='pop==0',
                    help="Choose a continent to add")
    pGroup.addParam('area2', params.StringParam, label='Country: ', default='All',
                    condition='pop==0 and area1!="All" and len(area1.split(",")) == 1',
                    help="Choose a country to add. Use the wizard to choose among the continent countries")
    pGroup.addParam('area3', params.StringParam, label='Country population: ', default='All',
                    condition='pop==0 and area1!="All" and len(area1.split(",")) == 1 and '
                              'area2!="All" and len(area2.split(",")) == 1',
                    help="Choose a country population to add. Use the wizard to choose among the populations in the "
                         "selected country")

    pGroup.addParam('ethnicity', params.StringParam, label='Ethnicity: ', condition='pop==1', default='All',
                    help="Choose an ethnicity to add")

    pGroup.addParam('addArea', params.LabelParam, label='Add defined population: ',
                    help='Add selected area or ethnicity to be included in the analysis')

    sGroup = form.addGroup('Populations summary')
    sGroup.addParam('inAreas', params.TextParam, width=70, default='',
                    label='Populations summary: ',
                    help='Summary of the areas or ethnicities to be included in the analysis')

  def _insertAllSteps(self):
    self._insertFunctionStep(self.coverageStep)
    self._insertFunctionStep(self.createOutputStep)

  def coverageStep(self):
    inEpiFile = self.getInputEpitopes()
    selPops = self.getSelectedPopulations()
    oFile = self.getCoverageOutputFile()

    fullPopStr = '","'.join(selPops)

    coveArgs = f'-p "{fullPopStr}" -c {self.getEnumText("mhc")} -f {inEpiFile} > {oFile} '

    iedbPlugin.runPopulationCoverage(self, coveArgs)

  def createOutputStep(self):
    coveDic = self.parseResults(self.getCoverageOutputFile())['average']

    outROIs = self.inputSequenceROIs.get().createCopy(self._getPath(), copyInfo=True, copyItems=True)
    outROIs._coverage = params.String(coveDic['coverage'])
    outROIs._averageHit = params.Float(coveDic['average_hit'])
    outROIs._pc90 = params.Float(coveDic['pc90'])
    outROIs._coverageFile = params.String(self.getCoverageOutputFile())

    if len(outROIs) > 0:
      self._defineOutputs(outputROIs=outROIs)

  ##################### UTILS #####################

  def getCoverageOutputFile(self):
    return os.path.abspath(self._getExtraPath('coverage_results.tsv'))

  def getInputEpitopes(self):
    inFile = self._getExtraPath('inputEpitopes.tsv')
    with open(inFile, 'w') as f:
      for roi in self.inputSequenceROIs.get():
        alleles = getattr(roi, '_alleles').get().replace('/', ',')
        f.write(f'{roi.getROISequence()}\t{alleles}\n')
    return inFile

  def parseResults(self, oFile):
    oDic = {}
    with open(oFile) as f:
      [f.readline() for i in range(2)]
      for line in f:
        if line.strip():
          sline = line.split('\t')
          oDic[sline[0]] = {'coverage': sline[1], 'average_hit': sline[2], 'pc90': sline[3].strip()}
        else:
          break
    return oDic

  def getAreaOptions(self, level):
    areas = [self.area1.get(), self.area2.get()]
    prevAreas = areas[:level-1]
    cDic = POP_DIC['Area']
    for pArea in prevAreas:
      cDic = cDic[pArea]

    pops = [area for area in cDic]
    pops.sort()
    return pops

  def getPopulationsElement(self):
    newArea = self.getEnumText('pop')
    if newArea == 'Area':
      areaLevel = 1
      cArea = getattr(self, f'area{areaLevel}').get()
      if cArea != 'All':
        newArea += f' / {cArea}'
      while cArea != 'All' and areaLevel < 3 and len(cArea.split(',')) == 1:
        areaLevel += 1
        cArea = getattr(self, f'area{areaLevel}').get()
        if cArea != 'All':
          newArea += f' / {cArea}'

    elif newArea == 'Ethnicity':
      if self.ethnicity.get() != 'All':
        newArea += f' / {self.ethnicity.get()}'
    return newArea

  def getSelectedPopulations(self):
    pops = []
    for popLine in self.inAreas.get().split('\n'):
      popElements = popLine.split('/ ')[-1].split(',')
      [pops.append(p.strip()) for p in popElements if p.strip() not in ['All', '']]

    if 'Area' in pops:
      pops.remove('Area')
      pops += list(POP_DIC['Area'].keys())

    pops.sort()
    return pops

  def _summary(self):
    summ = []
    if os.path.exists(self.getCoverageOutputFile()):
      coveDic = self.parseResults(self.getCoverageOutputFile())['average']
      summ.append(f'Average coverage of the set of sequence ROIs (MHC {self.getEnumText("mhc")} epitopes) for the selected populations: '
                  f'\n  - Coverage: {coveDic["coverage"]}'
                  f'\n  - Average hit: {coveDic["average_hit"]}'
                  f'\n  - PC90: {coveDic["pc90"]}'
                  f'\n\nFor more details, check the output file: {self.getCoverageOutputFile()}')
    return summ


