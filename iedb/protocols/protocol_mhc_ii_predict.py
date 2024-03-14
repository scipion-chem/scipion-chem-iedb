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
from ..constants import MHCII_alleles_dic
from ..utils import getAllMHCIIAlleles

class ProtMHCIIPrediction(EMProtocol):
  """Run a prediction using mhc-ii package from IEDB to predict MHC-II epitopes"""
  _label = 'mhc-ii prediction'

  _mhciMethodsDic = {'IEDB recommended': 'netmhciipan_el', 'Consensus-2.2': 'consensus',
                     'NetMHCIIpan-4.1': 'NetMHCIIpan', 'NN_align-1.0': 'nn_align', 'SMM_align-1.1': 'smm_align',
                     'Combinatorial Library-1.1': 'comblib', 'Sturniolo': 'sturniolo'}

  _species = ['Human', 'Mouse']
  _alleleGroups = ['7-allele method', 'Most frequent 26', 'Custom']
  _selTypes = ['Percentile rank', 'Score', 'Number of alleles', 'Top x%', 'Top x']

  def __init__(self, **kwargs):
    EMProtocol.__init__(self, **kwargs)

  def _defineParams(self, form):
    form.addSection(label='Input')
    iGroup = form.addGroup('Input')
    iGroup.addParam('inputSequence', params.PointerParam, pointerClass="Sequence",
                    label='Input protein sequence: ',
                    help="Protein sequence to perform the screening on")

    pGroup = form.addGroup('Parameters')
    pGroup.addParam('method', params.EnumParam, label='Prediction method: ',
                    default=0, choices=list(self._mhciMethodsDic.keys()),
                    help="Prediction method to use for the MHC-II binding prediction over the protein sequence.")
    pGroup.addParam('specie', params.EnumParam, label='Host species: ', expertLevel=params.LEVEL_ADVANCED,
                    default=0, choices=self._species,
                    help="Host specie  to predict the MHC-II epitopes on.")
    pGroup.addParam('lengths', params.StringParam, label='Peptide lengths: ', default='15',
                    help="Available lengths to include in the analysis.")

    pGroup.addParam('alleleGroup', params.EnumParam, label='Select allele groups: ', condition='specie==0',
                    default=1, choices=self._alleleGroups,
                    help="Select usual groups of human alleles to perform the analysis on")
    pGroup.addParam('alleleCustom', params.StringParam, label='Select alleles to add: ',
                    condition='specie!=0 or alleleGroup==2',
                    help="Host specie  to predict the MHC-II epitopes on.")

    sGroup = form.addGroup('Selection')
    sGroup.addParam('coreGroup', params.BooleanParam, label='Group peptides by core: ', default=True,
                    expertLevel=params.LEVEL_ADVANCED,
                    help="Whether to group the predicted peptides with the same 9mer core. The one with the better "
                         "score will be chosen as representant")
    sGroup.addParam('selType', params.EnumParam, label='Select output peptides by: ',
                    default=0, choices=self._selTypes,
                    help="Select output peptides in the chosen manner")
    sGroup.addParam('rank', params.FloatParam, label='Percentile rank threshold: ', default=10, condition='selType==0',
                    help="Predicted percentile rank threshold (<=)")
    sGroup.addParam('score', params.FloatParam, label='Score threshold: ', default=0.15, condition='selType==1',
                    help="Predicted percentile rank threshold (<=)")
    sGroup.addParam('topBinder', params.IntParam, label='Top threshold: ', default=5, condition='selType==2',
                    help="Select top x peptides based on percentile rank")
    sGroup.addParam('topPerc', params.FloatParam, label='Top percentage threshold (%): ', default=2,
                    condition='selType==3', help="Select top x% peptides based on percentile rank")
    sGroup.addParam('topN', params.IntParam, label='Top threshold: ', default=5, condition='selType==4',
                    help="Select top x peptides based on percentile rank")

    pGroup.addParam('remDup', params.BooleanParam, label='Remove duplicate peptides: ', default=True,
                    expertLevel=params.LEVEL_ADVANCED,
                    help="Whether to remove duplicate predicted peptides from the output")


  def _insertAllSteps(self):
    self._insertFunctionStep(self.mhcStep)
    self._insertFunctionStep(self.createOutputStep)

  def getAvailableAlleles(self):
    methKey = self._mhciMethodsDic[self.getEnumText('method')]
    alleDic = getAllMHCIIAlleles(methKey)
    return list(alleDic.keys())

  def writeInputFasta(self):
    faFile = self._getExtraPath('inputSequence.fa')
    self.inputSequence.get().exportToFile(faFile)
    return os.path.abspath(faFile)

  def getSelectedAlleles(self):
    '''Get list of alleles to perform the analysis on'''
    if self.specie.get() != 0 or self.alleleGroup.get() == 2:
      alList = self.alleleCustom.get().strip().split(', ')
    else:
      alList = MHCII_alleles_dic[self.getEnumText('alleleGroup')]
    return alList

  def filterAlleles(self, alDic, alList, lenList):
    '''Filters the allowed alleles and lengths from alList and lenList according to alDic and return the allele and
    length lists necessary to run predict_binding.py'''
    fAL, fLL = [], []
    for allele in alList:
      for length in lenList:
        if length in alDic[allele]:
          fAL.append(allele), fLL.append(length)
    return fAL, fLL

  def mhcStep(self):
    inFile = self.writeInputFasta()
    oFile = self.getMHCOutputFile()

    method = self._mhciMethodsDic[self.getEnumText('method')]
    selAlleles = self.getSelectedAlleles()
    # todo: lens must be between 11 and 30
    lenList = self.lengths.get().strip().split(', ')

    allowedAlleles = getAllMHCIIAlleles(method)
    alList = [allele for allele in selAlleles if allele in allowedAlleles]
    fullAlStr, fullLenStr = ','.join(alList), ','.join(lenList)

    mhcArgs = f'{method} {fullAlStr} {inFile} {fullLenStr} > {oFile} '

    iedbPlugin.runMHC_II(self, mhcArgs)

  def getCoreData(self, coreDic):
    cData = {'Epitope': [], 'Alleles': set([]), 'Score': 0 if self.selType == 1 else 100}
    for (idx, epitope) in coreDic:
      for (allele, score) in coreDic[(idx, epitope)]:
        cData['Alleles'].add(allele)
        if (score > cData['Score'] and self.selType == 1) or (score < cData['Score'] and self.selType != 1):
          cData['Score'], cData['Epitope'] = score, (idx, epitope)
    return cData

  def createOutputStep(self):
    epiDic = self.parseResults(self.getMHCOutputFile())

    inpSeq = self.inputSequence.get()
    outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
    method = self.getEnumText('method')

    for core in epiDic:
      if self.coreGroup.get():
        coreData = self.getCoreData(epiDic[core])
        epitope = coreData['Epitope'][1]
        idxs = [coreData['Epitope'][0], coreData['Epitope'][0] + len(epitope) - 1]
        roiSeq = Sequence(sequence=epitope, name='ROI_{}-{}'.format(*idxs), id='ROI_{}-{}'.format(*idxs),
                          description=f'MHC-II TepiTool epitope')
        seqROI = SequenceROI(sequence=inpSeq, seqROI=roiSeq, roiIdx=idxs[0], roiIdx2=idxs[1])
        seqROI._alleles, seqROI._core = params.String('/'.join(coreData['Alleles'])), params.String(core)
        seqROI._epitopeType = params.String('MHC-II')
        seqROI._source = params.String(method)
        setattr(seqROI, method, params.Float(coreData["Score"]))
        outROIs.append(seqROI)
      else:
        for (idxI, epitope) in epiDic[core]:
          idxs = [idxI, idxI + len(epitope) - 1]
          roiSeq = Sequence(sequence=epitope, name='ROI_{}-{}'.format(*idxs), id='ROI_{}-{}'.format(*idxs),
                            description=f'MHC-II TepiTool epitope')

          alleles, scores = [al[0] for al in epiDic[core][(idxI, epitope)]], \
                            [al[1] for al in epiDic[core][(idxI, epitope)]]
          allele, score = '/'.join(alleles), max(scores) if self.selType == 1 else min(scores)
          seqROI = SequenceROI(sequence=inpSeq, seqROI=roiSeq, roiIdx=idxs[0], roiIdx2=idxs[1])
          seqROI._alleles, seqROI._core = params.String(allele), params.String(core)
          seqROI._epitopeType = params.String('MHC-II')
          seqROI._source = params.String(self.method.get())
          setattr(seqROI, self.method.get(), params.Float(score))
          outROIs.append(seqROI)

    if len(outROIs) > 0:
      self._defineOutputs(outputROIs=outROIs)

  ##################### UTILS #####################

  def getMHCOutputFile(self):
    return os.path.abspath(self._getExtraPath('mhc-I_results.tsv'))

  def parseResults(self, oFile):
    '''Parse the results in the raw_output.tsv file generated by TepiTools and returns a dictionary
    as {allele: {(position, epitopeString): meanScore}}
    '''
    thType = self.selType.get()
    resAr = np.loadtxt(open(oFile, "rb"), delimiter="\t", skiprows=1, dtype=str)
    rankIdx = 8

    # Get the selection of the total results
    if thType == 0:
      scores = resAr[:, rankIdx].astype(float)
      idxSel = scores < self.rank.get()
    elif thType == 1:
      rankIdx = 7
      scores = resAr[:, rankIdx].astype(float)
      idxSel = scores > self.score.get()
    else:
      if thType == 2:
        nTop = int(self.topPerc.get() * 0.01 * resAr.shape[0])
      elif thType == 3:
        nTop = self.topN.get()
      idxSel = np.argsort(resAr, axis=0)
      idxSel = idxSel[:, rankIdx] < nTop

    resAr = resAr[idxSel, :]

    # Build the output from the selection
    epiDic = {}
    for row in resAr:
      allele, seq_id, pos, _, _, core, peptide, score, rank = row[:9]
      key = (int(pos), peptide)
      if not core in epiDic:
        epiDic[core] = {}
      if not key in epiDic[core]:
        epiDic[core][key] = []
      epiDic[core][key].append([allele, float(row[rankIdx])])
    return epiDic