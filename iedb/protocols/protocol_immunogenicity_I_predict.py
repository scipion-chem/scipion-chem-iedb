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

from pwem.protocols import EMProtocol
from pyworkflow.protocol import params

from pwchem.objects import SetOfSequenceROIs

from .. import Plugin as iedbPlugin

def writeInput(inSeqs, outFile):
  with open(outFile, 'w') as f:
    for seqId, seq in inSeqs.items():
      f.write(f'{seq}\n')
  return outFile

class ProtImmunogenicityPrediction(EMProtocol):
  """Run a prediction using immunogenicity package from IEDB to predict the immunogenicity of a peptide
  MHC (pMHC) complex over a set of sequence ROIs
  
  User IA Manual: ImmunogenicityIPredict Protocol

The ImmunogenicityIPredict protocol estimates the immunogenic potential of
MHC class I peptide binders. It is designed to help identify which peptide?MHC
complexes are more likely to elicit a T-cell response, going beyond simple
binding affinity predictions by incorporating immunogenicity-specific features.
This makes it particularly useful in the context of neoantigen discovery,
cancer immunotherapy, and vaccine design.

To begin, the user must provide a set of peptide sequences that are known or
predicted to bind to MHC class I molecules. These peptides should be in
FASTA format or imported through upstream protocols within the Scipion-Chem
platform. Each peptide is evaluated using a logistic regression model trained
on experimentally validated immunogenic and non-immunogenic MHC I ligands.

The protocol allows selection of the scoring method associated with the
immunogenicity model. The default method assigns a probability score between
zero and one to each peptide, reflecting the likelihood that the peptide will
trigger a cytotoxic T-cell response once presented by an MHC I molecule. Users
can set a classification threshold to distinguish between immunogenic and
non-immunogenic candidates, depending on the desired balance between sensitivity
and specificity.

Output from the protocol includes a prediction table in which each peptide is
assigned a numerical immunogenicity score, along with its classification status
based on the chosen threshold. These results can be ranked, filtered, or
visualized for comparison and prioritization. The protocol also supports
downstream integration with MHC binding, epitope clustering, or antigen selection
tools within the platform.

In summary, the ImmunogenicityIPredict protocol provides a statistically grounded,
sequence-based approach to evaluating MHC I immunogenicity. It supports the
identification of promising epitope candidates for immunotherapeutic and
vaccine-related applications, offering reproducible results and seamless
integration into computational immunology workflows.
  """
  _label = 'immunogenicity prediction'

  def __init__(self, **kwargs):
    EMProtocol.__init__(self, **kwargs)

  def getAvailableLengthList(self):
    return list(map(str, range(8, 15)))

  def _defineParams(self, form):
    form.addSection(label='Input')
    iGroup = form.addGroup('Input')
    iGroup.addParam('inputROIs', params.PointerParam, pointerClass="SetOfSequenceROIs",
                    label='Input sequence ROIs: ', 
                    help="Set of sequence ROIs to label with the present MHC-II alleles")

    iGroup.addParam('customMask', params.StringParam, label='Anchor positions: ', default='',
                    expertLevel=params.LEVEL_ADVANCED,
                    help="A custom mask scheme can be entered by inputting comma separated numbers."
                         "As a default, the first, second, and C-terminal amino acid are masked")

  def _insertAllSteps(self):
    self._insertFunctionStep(self.immunoStep)
    self._insertFunctionStep(self.createOutputStep)
  
  def getInputSequences(self):
    '''Returns the input sequences as a dict like {seqId: sequence}
    '''
    inSeqs = {}
    inSet = self.inputROIs.get()
    for item in inSet:
      seq, name = item.getROISequence(), item.getROIId()
      inSeqs[name] = seq
    return inSeqs

  
  def immunoStep(self):
    inSeqs = self.getInputSequences()
    inFile = writeInput(inSeqs, self._getExtraPath('inputSequences.txt'))
    oFile = self.getImmunoOutputFile()

    immunoArgs = f'{inFile} > {oFile} '
    iedbPlugin.runImmunogenicity(self, immunoArgs)

  def createOutputStep(self):
    epiDic = self.parseResults(self.getImmunoOutputFile())
    outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))

    # Each input ROI is labelled with the alleles found inside them
    for curROI in self.inputROIs.get():
      roiSeq = curROI.getROISequence()
      curROI._immunogenicity = params.Float(epiDic[roiSeq])
      outROIs.append(curROI)

    if len(outROIs) > 0:
      self._defineOutputs(outputROIs=outROIs)

  ##################### UTILS #####################

  def getImmunoOutputFile(self):
    return os.path.abspath(self._getExtraPath('immunogenicity_results.tsv'))

  def parseResults(self, oFile):
    '''Parse the results in the output file generated by Immunogenicity and returns a dictionary
    as {sequence: score}
    '''
    # Build the output from the selection
    epiDic, start = {}, False
    with open(oFile) as f:
      for line in f:
        if start:
          sline = line.split(',')
          epiDic[sline[0]] = float(sline[2])
        
        if line.strip() == 'peptide,length,score':
          start = True
        
    return epiDic