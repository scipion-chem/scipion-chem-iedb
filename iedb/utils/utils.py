# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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

from iedb import Plugin
from iedb.constants import POP_DIC

def getAllelesFile(mhc, method):
  return Plugin.getPluginHome(f'constants/alleles-{mhc}/{method}_alleles.txt')

def getAllMHCIAlleles(method, specie='human'):
  '''Parse the possible alleles given a method from the files stored in constants.
  Returns a dictionary of the form: {specie: {mhc: [lengths]}}'''
  alDic = {}
  alleFile = getAllelesFile('I', method)
  specie = specie.lower()
  with open(alleFile) as f:
    f.readlines(2)
    for line in f:
      sp, allele, l = line.split()
      if specie == sp:
        if allele not in alDic:
          alDic[allele] = []

        alDic[allele].append(l)
  return alDic

def getAllMHCIIAlleles(method, specie='human'):
  alleles = []
  if specie == 'human':
    alleFile = getAllelesFile('II', method)
    with open(alleFile) as f:
      f.readline()
      for line in f:
        alleles.append(line.strip())
  else:
    from ..constants import MOUSE_MHCII_ALLELES
    alleles = MOUSE_MHCII_ALLELES
  return alleles

def getMHCAlleles(roi, mhc):
  alleles = []
  if hasattr(roi, '_allelesMHCI') and mhc in ['I', 'combined']:
    alleles += [getattr(roi, '_allelesMHCI').get().replace('/', ',')]

  if hasattr(roi, '_allelesMHCII') and mhc in ['II', 'combined']:
    alleles += [getattr(roi, '_allelesMHCII').get().replace('/', ',')]
  return ','.join(alleles)

def writeInputEpitopeFiles(inputROIs, epFile, separated, mhc):
  outFiles = []
  if separated:
    for roi in inputROIs:
      groupFile = epFile.replace('.tsv', f'_{roi.clone().getObjId()}.tsv')
      with open(groupFile, 'w') as f:
        alleles = getMHCAlleles(roi, mhc)
        f.write(f'{roi.getROISequence()}\t{alleles}\n')
      outFiles.append(groupFile)
  else:
    epFile = epFile.replace('.tsv', f'_All.tsv')
    with open(epFile, 'w') as f:
      for roi in inputROIs:
        alleles = getMHCAlleles(roi, mhc)
        f.write(f'{roi.getROISequence()}\t{alleles}\n')
    outFiles.append(epFile)

  return outFiles

def translateArea(pops):
  if 'Area' in pops:
    pops.remove('Area')
    pops += list(POP_DIC['Area'].keys())

  pops.sort()
  return pops

def buildMHCCoverageArgs(inputROIs, epFile, populations, mhc, oDir, separated=True):
  inEpiFiles = writeInputEpitopeFiles(inputROIs, epFile, separated, mhc)
  fullPopStr = '","'.join(populations)

  coveArgs = []
  for epFile in inEpiFiles:
    epBase = os.path.basename(epFile)
    oFile = os.path.join(oDir, epBase.replace('.tsv', '_results.tsv'))
    coveArgs += [f'-p "{fullPopStr}" -c {mhc} -f {epFile} > {oFile} ']
  return coveArgs

def parseCoverageResults(oFile):
  oDic = {}
  with open(oFile) as f:
    [f.readline() for i in range(2)]
    for line in f:
      if line.strip():
        sline = line.split('\t')
        oDic[sline[0]] = {'coverage': float(sline[1][:-2]), 'average_hit': sline[2], 'pc90': sline[3].strip()}
      else:
        break
  return oDic