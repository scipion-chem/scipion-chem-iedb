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

SEQ, SEQROIS = 0, 1
RANK, SCORE, TOPP, NTOP = 0, 1, 2, 3

class ProtMHCIIPrediction(EMProtocol):
  """Run a prediction using mhc-ii package from IEDB to predict MHC-II epitopes
  
  User Manual: MhcIIPredict Protocol in Scipion-Chem-IEDB

The MhcIIPredict protocol allows users to predict the binding affinity of peptides
to MHC class II molecules using models available through the IEDB resource. Unlike
MHC class I, where peptides typically bind with fixed length, MHC class II molecules
can accommodate longer peptides with variable core binding regions. This protocol
automates the identification and scoring of such core sequences within user-supplied
peptides.

To begin, the user must provide a set of peptides, usually in FASTA format or as
output from a previous protocol. These peptides are typically between 13 and 25
amino acids long. The protocol automatically scans overlapping 9-mer cores within
each peptide and evaluates their predicted affinity to selected MHC class II alleles.

The user may choose from a broad range of MHC II alleles, including human HLA-DR,
DP, and DQ types, as well as murine variants. Multiple alleles can be selected
simultaneously to evaluate peptide promiscuity or allele-specific binding. The
prediction method must also be selected. Available options include recommended
tools such as NetMHCIIpan and consensus-based methods supported by the IEDB
infrastructure.

Each peptide is assigned a quantitative binding score, which may be expressed
as predicted affinity or percentile rank, depending on the selected method. A
threshold can be configured to classify peptides as binders or non-binders, with
additional labels for strong or weak binders. These labels assist in prioritizing
epitopes for vaccine development or immunogenicity testing.

The output consists of a table summarizing each peptide?allele pair with its
best scoring core, associated binding score, and classification. The protocol
can optionally highlight the predicted binding core within the full-length
peptide, which is useful for experimental validation or peptide synthesis
design. All results are compatible with downstream Scipion-Chem protocols,
including those for immunogenicity scoring, epitope clustering, or structural
mapping.

In summary, the MhcIIPredict protocol offers a robust and flexible interface for
predicting peptide?MHC II interactions. It supports multiple alleles and
prediction methods, and it enables efficient screening of long peptides for
potential T-helper cell epitopes within a reproducible and integrated
bioinformatics environment.
  
  """
  _label = 'mhc-ii prediction'

  MINLEN, MAXLEN = 11, 30
  selMap = {RANK: 'rank', SCORE: 'score', TOPP: 'topPerc', NTOP: 'topN'}

  _mhciiMethodsDic = {'IEDB recommended': 'netmhciipan_el', 'Consensus-2.2': 'consensus',
                     'NetMHCIIpan-4.1': 'NetMHCIIpan', 'NN_align-1.0': 'nn_align', 'SMM_align-1.1': 'smm_align',
                     'Combinatorial Library-1.1': 'comblib', 'Sturniolo': 'sturniolo'}
  _species = ['Human', 'Mouse']
  _alleleGroups = ['7-allele method', 'Most frequent 26', 'Custom']
  _selTypes = ['Percentile rank', 'Score', 'Top x%', 'Top x']

  def __init__(self, **kwargs):
    super().__init__(**kwargs)

  def getAvailableLengthList(self):
    return list(map(str, range(self.MINLEN, self.MAXLEN)))
  
  def _defineInputParams(self, form):
    iGroup = form.addGroup('Input')
    iGroup.addParam('inputSource', params.EnumParam, label='Input action: ',
                    default=SEQ, choices=['Predict over sequence', 'Label sequence ROIs'],
                    display=params.EnumParam.DISPLAY_HLIST,
                    help="Whether to predict the MHC-I epitopes present in a sequence and label them with the present"
                         "alleles or label an already existing SetOfSequenceROIs with the alleles found in them")
    iGroup.addParam('inputSequence', params.PointerParam, pointerClass="Sequence", allowsNull=True,
                    label='Input protein sequence: ', condition='inputSource==0',
                    help="Protein sequence to perform the screening on")
    iGroup.addParam('inputSequenceROIs', params.PointerParam, pointerClass="SetOfSequenceROIs",
                    label='Input sequence ROIs: ', condition='inputSource==1', allowsNull=True,
                    help="Set of sequence ROIs to label with the present MHC-II alleles")
    return iGroup

  def _defineParams(self, form):
    form.addSection(label='Input')
    iGroup = self._defineInputParams(form)
    
    pGroup = form.addGroup('Parameters')
    pGroup.addParam('method', params.EnumParam, label='Prediction method: ',
                    default=0, choices=list(self._mhciiMethodsDic.keys()),
                    help="Prediction method to use for the MHC-II binding prediction over the protein sequence.")
    pGroup.addParam('specie', params.EnumParam, label='Host species: ', expertLevel=params.LEVEL_ADVANCED,
                    default=0, choices=self._species,
                    help="Host specie  to predict the MHC-II epitopes on.")
    pGroup.addParam('lengths', params.StringParam, label='Peptide lengths: ', default='15',
                    help="Available lengths to include in the analysis. "
                         "You can include several lengths as comma separated  (11, 12,13)"
                         f"Lenght limits are [{self.MINLEN}, {self.MAXLEN}]")

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
    sGroup.addParam('mergeAlleles', params.BooleanParam, label='Merge alleles: ', default=True,
                    expertLevel=params.LEVEL_ADVANCED, condition=f'inputSource=={SEQ}',
                    help="Merges same epitope sequences predicted to interact with different alleles into the same "
                         "sequence ROI")

    sGroup.addParam('selType', params.EnumParam, label='Select output peptides by: ',
                    default=0, choices=self._selTypes,
                    help="Select output peptides in the chosen manner")
    sGroup.addParam('rank', params.FloatParam, label='Percentile rank threshold: ', default=10, condition='selType==0',
                    help="Predicted percentile rank threshold (<=)")
    sGroup.addParam('score', params.FloatParam, label='Score threshold: ', default=0.15, condition='selType==1',
                    help="Predicted percentile rank threshold (<=)")
    sGroup.addParam('topPerc', params.FloatParam, label='Top percentage threshold (%): ', default=2,
                    condition='selType==2', help="Select top x% peptides based on percentile rank")
    sGroup.addParam('topN', params.IntParam, label='Top threshold: ', default=5, condition='selType==3',
                    help="Select top x peptides based on percentile rank")



  def _insertAllSteps(self):
    self._insertFunctionStep(self.mhcStep)
    self._insertFunctionStep(self.createOutputStep)

  def getLenghts(self):
    return [int(leni.strip()) for leni in self.lengths.get().strip().split(',')]

  def mhcStep(self):
    inFile = self.writeInputFasta()
    oFile = self.getMHCOutputFile()

    method = self._mhciiMethodsDic[self.getEnumText('method')]
    selAlleles = self.getSelectedAlleles()
    lenList = self.getLenghts()

    allowedAlleles = getAllMHCIIAlleles(method)
    alList = [allele for allele in selAlleles if allele in allowedAlleles]
    fullAlStr, fullLenStr = ','.join(alList), ','.join([str(l) for l in lenList])

    mhcArgs = f'{method} {fullAlStr} {inFile} {fullLenStr} > {oFile} '

    iedbPlugin.runMHC_II(self, mhcArgs)

  def createOutputStep(self):
    epiDic = self.parseResults(self.getMHCOutputFile())

    inpSeq = self.inputSequence.get()
    outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
    method = f"MHCII_{self.getEnumText('method')}"

    if self.inputSource.get() == SEQ:
      epiDic, epitopesList = epiDic['1'], []
      for core in epiDic:
        epitopesList += self.mergeCoreData(epiDic[core], merge=self.coreGroup.get())

      if self.mergeAlleles.get():
        epitopesList = self.mergeAllelesEpitopes(epitopesList)

      for (idxI, epitope), alleles, score in epitopesList:
        idxs = [idxI, idxI + len(epitope) - 1]
        roiSeq = Sequence(sequence=epitope, name='ROI_{}-{}'.format(*idxs), id='ROI_{}-{}'.format(*idxs),
                          description=f'MHC-II TepiTool epitope')

        seqROI = SequenceROI(sequence=inpSeq, seqROI=roiSeq, roiIdx=idxs[0], roiIdx2=idxs[1])
        seqROI._allelesMHCII = params.String('/'.join(alleles))
        seqROI._epitopeType = params.String('MHC-II')
        seqROI._source = params.String(method)
        setattr(seqROI, method, params.Float(score))
        outROIs.append(seqROI)

    else:
      # Each input ROI is labelled with the alleles found inside them
      inROIs = [roi.clone() for roi in self.inputSequenceROIs.get()]
      i, lens = 0, self.getLenghts()
      for curROI in inROIs:
        curAlleles, curScores = [], []
        if len(curROI.getROISequence()) >= min(lens):
          i += 1
          roiId = str(i)
          if roiId in epiDic:
            for core in epiDic[roiId]:
              for (idxI, epitope) in epiDic[roiId][core]:
                alleles, scores = [al[0] for al in epiDic[roiId][core][(idxI, epitope)]], \
                                  [al[1] for al in epiDic[roiId][core][(idxI, epitope)]]
                curAlleles += alleles
                curScores += scores

        allele, score = '/'.join(curAlleles), min(curScores) if curScores else 0
        curROI._allelesMHCII = params.String(allele)
        curROI._sourceMHCII = params.String(method)
        setattr(curROI, method, params.Float(score))
        outROIs.append(curROI)

    if len(outROIs) > 0:
      self._defineOutputs(outputROIs=outROIs)

  ##################### UTILS #####################

  def getAvailableAlleles(self):
    methKey = self._mhciiMethodsDic[self.getEnumText('method')]
    alleDic = getAllMHCIIAlleles(methKey)
    return list(alleDic.keys())

  def writeInputFasta(self):
    faFile = self._getExtraPath('inputSequence.fa')
    if self.inputSource.get() == 0:
      self.inputSequence.get().exportToFile(faFile)
    else:
      lens = self.getLenghts()
      self.inputSequenceROIs.get().exportToFile(faFile, mainSeq=False, minLen=min(lens))
    return os.path.abspath(faFile)

  def getSelectedAlleles(self):
    '''Get list of alleles to perform the analysis on'''
    if self.specie.get() != 0 or self.alleleGroup.get() == 2:
      alList = self.alleleCustom.get().strip().split(', ')
    else:
      alList = MHCII_alleles_dic[self.getEnumText('alleleGroup')]
    return alList

  def mergeCoreData(self, coreDic, merge=True):
    '''Merges the epitopes containing the best core. Alleles are appended and the best score is taken.
    :param coreDic: dic containing the same core epitopes as {(idx, epitopeStr): [(alleles, score)]}
    :return: list of epitopes described as [ ((idx, epitopeStr), alleles, score), ... ]
    '''
    allEpitopes = []
    key, alleles, score = (None, None), set([]), 0 if self.selType == 1 else 100
    for newKey in coreDic:
      for (newAlleles, newScore) in coreDic[newKey]:
        if merge:
          alleles.add(newAlleles)
          if (newScore > score and self.selType == 1) or (newScore < score and self.selType != 1):
            score, key = newScore, newKey
        else:
          allEpitopes.append((newKey, newAlleles, newScore))
    if merge:
      allEpitopes.append((key, list(alleles), score))
    return allEpitopes

  def mergeAllelesEpitopes(self, epitopeList):
    '''Merge the data of those epitopes whose sequence is equal. Alleles are appended, best score is kept
    :param epitopeList: list of epitopes as [ ((idx, epitopeStr), alleles, score), ... ]
    :return: list of non-duplicated epitopes as [ ((idx, epitopeStr), alleles, score), ... ]
    '''
    uniEpitopes = {}
    for (key, newAlleles, newScore) in epitopeList:
       if key in uniEpitopes:
         alleles, score = uniEpitopes[key]
         newScore = newScore if (newScore > score and self.selType == 1) or (newScore < score and self.selType != 1) \
           else score
         newAlleles += alleles

       uniEpitopes[key] = (newAlleles, newScore)

    outEps = []
    for key, (alleles, score) in uniEpitopes.items():
      alleles = list(set(alleles))
      outEps.append((key, alleles, score))

    return outEps

  def getMHCOutputFile(self):
    return os.path.abspath(self._getExtraPath('mhc-II_results.tsv'))

  def getRankIdx(self):
    return 7 if self.selType.get() == SCORE else 8

  def getResultsArray(self, oFile):
    resAr = np.loadtxt(open(oFile, "rb"), delimiter="\t", skiprows=1, dtype=str)

    thType = self.selType.get()
    rankIdx = self.getRankIdx()
    selValue = getattr(self, self.selMap[thType]).get()

    if thType in [RANK, SCORE]:
      scores = resAr[:, rankIdx].astype(float)
      idxSel = scores < selValue
    else:
      if thType == TOPP:
        nTop = int(selValue * 0.01 * resAr.shape[0])
      elif thType == NTOP:
        nTop = selValue
      idxSel = np.argsort(resAr, axis=0)
      idxSel = idxSel[:, rankIdx] < nTop

    return resAr[idxSel, :]

  def parseResults(self, oFile):
    '''Parse the results in the raw_output.tsv file generated by TepiTools and returns a dictionary
    as {seq_id: {core: {(position, epitopeString): [allele, score]}}}
    '''
    resAr = self.getResultsArray(oFile)

    # Build the output from the selection
    epiDic = {}
    for row in resAr:
      allele, seq_id, pos, _, _, core, peptide, score, rank = row[:9]
      key = (int(pos), peptide)
      if not seq_id in epiDic:
        epiDic[seq_id] = {}
      if not core in epiDic[seq_id]:
        epiDic[seq_id][core] = {}
      if not key in epiDic[seq_id][core]:
        epiDic[seq_id][core][key] = []

      rankIdx = self.getRankIdx()
      epiDic[seq_id][core][key].append([allele, float(row[rankIdx])])
    return epiDic

  def getInputSequences(self):
    if self.inputSource.get() == SEQ:
      return [self.inputSequence.get().getSequence()]
    else:
      return [roi.getROISequence() for roi in self.inputSequenceROIs.get()]

  def _validate(self):
    vals = []
    lens = self.getLenghts()
    if min(lens) < self.MINLEN or max(lens) > self.MAXLEN:
      vals.append(f'Length of the epitopes must be between {self.MINLEN} and {self.MAXLEN} both included')

    inSeqs = self.getInputSequences()
    if self.inputSource.get() == SEQ:
      for s in inSeqs:
        if len(s) < min(lens):
          vals.append(f'Input sequences must be at least {min(lens)} aminoacids long '
                      f'(The smallest length you have defined). Please check your input.')
          break
    elif self.inputSource.get() == SEQROIS:
      c = 0
      for s in inSeqs:
        if len(s) >= min(lens):
          c += 1

      if c == 0:
        vals.append(f'None of the input sequence ROIs is at least {min(lens)} residues long, so the analysis cannot '
                    f'be performed. Please check your input.')

    return vals

  def _warnings(self):
    warns = []
    if self.inputSource.get() == SEQROIS:
      lens = self.getLenghts()
      inSeqs = self.getInputSequences()
      for s in inSeqs:
        if len(s) < min(lens):
          warns.append(f'At least one of your sequences is not {min(lens)} aminoacids long '
                       f'(The smallest length you have defined). These sequences will have Null scores and alleles')
          break
    return warns
