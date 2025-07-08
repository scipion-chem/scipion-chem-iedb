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

import os, json
from Bio.PDB.PDBParser import PDBParser

from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pwem.convert import AtomicStructHandler

from pwchem.objects import Sequence, SequenceROI, SetOfSequenceROIs, StructROI, SetOfStructROIs
from pwchem.utils import runOpenBabel, createPocketFile

from .. import Plugin as iedbPlugin
from ..constants import ELLI_DIC

class ProtElliProPrediction(EMProtocol):
  """Run a prediction using ElliPro to extract B-cell structural epitopes
  
  User IA Manual: ElliproPredict Protocol

The ElliproPredict protocol is designed to predict discontinuous B-cell epitopes
based on the three-dimensional structure of a protein. It uses the ElliPro method,
which combines residue protrusion indices with clustering algorithms to identify
surface-exposed regions likely to be recognized by antibodies. The protocol is
particularly useful for mapping conformational epitopes that cannot be inferred
from sequence alone.

To run the protocol, the user must provide a protein structure file in PDB format.
The structure should be complete in the regions of interest, as missing residues
can affect the accuracy of the surface and shape-based predictions. The algorithm
requires a well-defined 3D model to compute the geometric features associated
with potential epitopic sites.

The user can adjust the minimum protrusion index threshold, which defines how far
a residue must extend from the molecular surface to be considered as part of an
epitope. A higher threshold will result in fewer, more conservative predictions.
Additionally, the user may define the clustering distance to control how residues
are grouped into predicted epitopic regions. This parameter affects the size and
compactness of the resulting epitope candidates.

The protocol outputs a list of predicted epitope clusters, each associated with a
set of residues and a confidence score derived from the average protrusion values.
These results are presented in a structured format and can be visualized on the
protein structure, allowing users to inspect the spatial distribution and
accessibility of the predicted sites.

In summary, the ElliproPredict protocol enables structure-based identification of
conformational B-cell epitopes, providing valuable insight into antigenic regions
that may not be evident from sequence-based methods. It supports applications in
vaccine design, antibody engineering, and immunogenicity assessment, offering a
reproducible and interpretable tool within the Scipion-Chem platform.
  """
  _label = 'ellipro prediction'

  def __init__(self, **kwargs):
    EMProtocol.__init__(self, **kwargs)

  def _defineParams(self, form):
    form.addSection(label='Input')
    iGroup = form.addGroup('Input')
    iGroup.addParam('inputAtomStruct', params.PointerParam, pointerClass="AtomStruct",
                    label='Input protein structure: ',
                    help="Protein structure to perform the prediction on")

    iGroup.addParam("rchains", params.BooleanParam,
                    label='Select specific chains: ', default=False, important=True,
                    help='Keep only the chains selected')
    iGroup.addParam("chain_name", params.StringParam,
                    label="Keep chains: ", important=True, condition="rchains",
                    help="Select the chain(s) you want to use in the analysis")

    pGroup = form.addGroup('Parameters')
    pGroup.addParam('minScore', params.FloatParam, label='Minimum score: ', default=0.5,
                    help="Minimum score for a predicted epitope to be considered")
    pGroup.addParam('maxDist', params.FloatParam, label='Maximum distance (A): ', default=6.0,
                    help="Maximum distance between the epitope points")


  def _insertAllSteps(self):
    self._insertFunctionStep(self.convertInputStep)
    self._insertFunctionStep(self.elliProStep)
    self._insertFunctionStep(self.createOutputStep)

  def convertInputStep(self):
    inASFile = self.inputAtomStruct.get().getFileName()
    if inASFile.endswith('.pdb'):
      os.link(inASFile, self.getReceptorPDB())
    else:
      self.convertReceptor2PDB(inASFile)

  def elliProStep(self):
    chains = self.getInputChains()
    elliArgs = f' -f {self.getReceptorPDB()} -s {self.minScore.get()} -d {self.maxDist.get()} '
    if chains:
      elliArgs += f'-c {",".join(chains)}'
    iedbPlugin.runElliPro(self, elliArgs, cwd=os.path.abspath(self._getExtraPath()))

  def createOutputStep(self):
    epiDic = self.parseResults()
    handler = AtomicStructHandler(self.getReceptorPDB())
    structModel = PDBParser().get_structure('inputAS', self.getReceptorPDB())[0]

    outSeqROIOuts = {}
    for chainId, chainDic in epiDic['Linear'].items():
      mapDic = self.getResidueIdMap(structModel, chainId)

      outSeqROIOuts[chainId] = SetOfSequenceROIs(filename=self._getPath(f'sequenceROIs_{chainId}.sqlite'))
      chainSeq = str(handler.getSequenceFromChain(modelID=0, chainID=chainId))
      outSeq = Sequence(name=f'Chain {chainId}', sequence=chainSeq, description=f'Chain {chainId}')

      for roiId, roiDic in chainDic.items():
        idxs = roiDic['Index']
        idxs = [mapDic[i] for i in idxs]
        roiSeq = Sequence(sequence=roiDic['Peptide'], name='ROI_{}-{}'.format(*idxs), id='ROI_{}-{}'.format(*idxs),
                          description='Linear ElliPro epitope')
        seqROI = SequenceROI(sequence=outSeq, seqROI=roiSeq, roiIdx=idxs[0], roiIdx2=idxs[1])
        seqROI._epitopeType = params.String('B')
        seqROI._source = params.String('ElliPro')
        seqROI._type = 'ElliPro'
        setattr(seqROI, 'ElliPro', params.Float(roiDic['Score']))
        outSeqROIOuts[chainId].append(seqROI)

    for chainId, outROIs in outSeqROIOuts.items():
      self._defineOutputs(**{f'outSequenceROIs_{chainId}': outROIs})


    outROIs = SetOfStructROIs(filename=self._getPath('StructROIs.sqlite'))
    for roiId, roiDic in epiDic['Discontinous'].items():
      roiCoords = []
      for res in roiDic['Residues']:
        chainId, resIdx = res.split(':')[0], int(res.split(':')[1][1:])
        residue = structModel[chainId][resIdx]
        atoms = residue.get_atoms()
        for a in atoms:
          roiCoords.append(list(a.get_coord()))

      pocketFile = createPocketFile(roiCoords, roiId, outDir=self._getExtraPath())
      pocket = StructROI(pocketFile, self.getReceptorPDB())
      pocket.calculateContacts()
      outROIs.append(pocket)

    if len(outROIs) > 0:
      self._defineOutputs(outputStructROIs=outROIs)

  ##################### UTILS #####################
  def getResidueIdMap(self, structModel, chainId):
    resIdxs = {}
    for i, res in enumerate(structModel[chainId]):
      resIdxs[res.get_id()[1]] = i+1
    return resIdxs

  def getReceptorPDB(self):
    return os.path.abspath(self._getExtraPath('inputAS.pdb'))

  def convertReceptor2PDB(self, proteinFile):
    inExt = os.path.splitext(os.path.basename(proteinFile))[1]
    oFile = self.getReceptorPDB()
    if not os.path.exists(oFile):
      args = ' -i{} {} -opdb -O {}'.format(inExt[1:], os.path.abspath(proteinFile), oFile)
      runOpenBabel(protocol=self, args=args, cwd=self._getTmpPath())
    return oFile

  def getInputChains(self):
    chains = []
    if self.rchains.get() and self.chain_name.get().strip():
      chDic = json.loads(self.chain_name.get())
      if 'residues' in chDic:
        chains = [chDic['chain'].strip()]
      elif 'model-chain' in chDic:
        mcs = chDic['model-chain'].split(',')
        chains = [mc.split('-')[1].strip() for mc in mcs]
    return chains

  def parseResults(self):
    rDic = {'Linear': {}, 'Discontinous': {}}
    with open(self._getExtraPath('output.txt')) as f:
      for line in f:
        line = line.strip()
        if line:
          sl = line.split(',')
          mode = sl[-1]
          if mode == 'Linear':
            chain = sl[2]
            if chain not in rDic['Linear']:
              rDic['Linear'][chain] = {}
            rDic['Linear'][chain][sl[0]] = {'Index': [int(sl[3]), int(sl[4])], 'Peptide': sl[5], 'Score': float(sl[7])}
          elif mode == 'Discontinous':
            rDic['Discontinous'][sl[0]] = {'Residues': [res.strip() for res in sl[2:-3]], 'Score': float(sl[-2])}

    return rDic



