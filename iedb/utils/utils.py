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

from iedb import Plugin

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