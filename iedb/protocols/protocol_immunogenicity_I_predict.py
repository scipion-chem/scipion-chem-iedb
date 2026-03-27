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
  
   AI Generated:

      ProtImmunogenicityPrediction - User Manual

      Overview
      --------
      The ProtImmunogenicityPrediction protocol evaluates the immunogenic
      potential of peptide sequences presented by MHC class I molecules. It
      predicts the likelihood that a peptide–MHC complex will elicit a
      cytotoxic T-cell response.

      This protocol goes beyond binding affinity prediction by incorporating
      sequence-derived features associated with T-cell recognition, making it
      particularly useful for neoantigen discovery, cancer immunotherapy, and
      vaccine design.

      Input Requirements
      ------------------
      1. **Peptide Sequences**:
         - A SetOfSequenceROIs containing peptide regions of interest.
         - Each ROI represents a candidate epitope sequence.

      2. **Sequence Characteristics**:
         - Typically peptides of length 8–14 amino acids.
         - Should correspond to predicted or known MHC class I binders.

      Workflow
      --------
      1. **Input Extraction**:
         - Peptide sequences are extracted from the input ROIs.
         - Sequences are written into a plain text input file.

      2. **Immunogenicity Prediction**:
         - Runs the IEDB immunogenicity prediction tool.
         - Applies a statistical model (logistic regression-based).
         - Computes a score representing immunogenic potential.

      3. **Scoring**:
         - Each peptide is assigned a numerical score.
         - Higher scores indicate higher likelihood of T-cell activation.

      4. **Annotation**:
         - Scores are mapped back to the original ROIs.
         - Each ROI is labeled with its immunogenicity value.

      Outputs
      -------
      - **Annotated Sequence ROIs (SetOfSequenceROIs)**:
        - Input ROIs enriched with immunogenicity scores.
        - Each ROI includes:
          - Original peptide sequence
          - Immunogenicity score
          - Metadata preserved from input

      - **Prediction Output File**:
        - Tabular file containing peptide sequences and scores.
        - Stored for traceability and further analysis.

      Advanced Options
      ----------------
      - Custom masking of anchor residues:
        - Allows specification of positions to ignore in scoring.
        - Useful for fine-tuning predictions based on MHC binding motifs.

      Validation & Warnings
      ---------------------
      - Input peptides should be valid MHC class I binders.
      - Very short or very long peptides may not be properly evaluated.
      - Duplicate sequences may overwrite scores if not handled carefully.
      - Results depend on the statistical model and training dataset.

      Practical Recommendations
      -------------------------
      - Use this protocol after MHC binding prediction steps.
      - Focus on peptides with strong binding affinity first.
      - Combine immunogenicity scores with other filters (e.g., expression data).
      - Validate top candidates experimentally when possible.

      Final Perspective
      -----------------
      ProtImmunogenicityPrediction provides a fast and effective method to
      prioritize peptide candidates based on their potential to trigger an
      immune response. By integrating immunogenicity prediction into the
      workflow, it enhances the selection of biologically relevant epitopes
      for therapeutic and vaccine applications.
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