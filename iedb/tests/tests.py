# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307 USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.tests import setupTestProject, DataSet, BaseTest

from pwem.protocols import ProtImportSequence, ProtImportPdb

from pwchem.utils import assertHandle

from ..protocols import *


class BaseImportSeq(BaseTest):
	NAME = 'USER_SEQ'
	DESCRIPTION = 'User description'
	AMINOACIDSSEQ1 = 'MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG'

	@classmethod
	def setUpClass(cls):
		super().setUpClass()
		cls.ds = DataSet.getDataSet('model_building_tutorial')
		setupTestProject(cls)

		cls._runImportSeq()
		cls._waitOutput(cls.protImportSeq, 'outputSequences', sleepTime=5)

	@classmethod
	def _runImportSeq(cls):
		kwargs = {'inputSequenceName': cls.NAME,
							'inputSequenceDescription': cls.DESCRIPTION,
							'inputRawSequence': cls.AMINOACIDSSEQ1
							}

		cls.protImportSeq = cls.newProtocol(
			ProtImportSequence, **kwargs)
		cls.proj.launchProtocol(cls.protImportSeq, wait=False)


class TestBepiPredPrediction(BaseImportSeq):
	def _runBepiPredPrediction(self):
		protBepiPred = self.newProtocol(ProtBepiPredPrediction)

		protBepiPred.inputSequence.set(self.protImportSeq)
		protBepiPred.inputSequence.setExtended('outputSequence')

		self.proj.launchProtocol(protBepiPred, wait=False)
		return protBepiPred

	def test(self):
		protBepiPred = self._runBepiPredPrediction()
		self._waitOutput(protBepiPred, 'outputROIs', sleepTime=10)
		assertHandle(self.assertIsNotNone, getattr(protBepiPred, 'outputROIs', None))

class TestMHCIPrediction(BaseImportSeq):
	def _runMHCIPrediction(self):
		protMHCI = self.newProtocol(ProtMHCIPrediction)

		protMHCI.inputSequence.set(self.protImportSeq)
		protMHCI.inputSequence.setExtended('outputSequence')

		self.proj.launchProtocol(protMHCI, wait=False)
		return protMHCI

	def test(self):
		protMHCI = self._runMHCIPrediction()
		self._waitOutput(protMHCI, 'outputROIs', sleepTime=10)
		assertHandle(self.assertIsNotNone, getattr(protMHCI, 'outputROIs', None))

class TestMHCIIPrediction(BaseImportSeq):
	def _runMHCIIPrediction(self):
		protMHCII = self.newProtocol(ProtMHCIIPrediction)

		protMHCII.inputSequence.set(self.protImportSeq)
		protMHCII.inputSequence.setExtended('outputSequence')

		self.proj.launchProtocol(protMHCII, wait=False)
		return protMHCII

	def test(self):
		protMHCII = self._runMHCIIPrediction()
		self._waitOutput(protMHCII, 'outputROIs', sleepTime=10)
		assertHandle(self.assertIsNotNone, getattr(protMHCII, 'outputROIs', None))


class TestMHCPopulationCoverage(TestMHCIIPrediction):
	def _runMHCCoverage(self, protMHC):
		protPop = self.newProtocol(ProtMHCIIPopulationCoverage,
															 mhc=1, inAreas='Area')

		protPop.inputSequenceROIs.set(protMHC)
		protPop.inputSequenceROIs.setExtended('outputROIs')

		self.proj.launchProtocol(protPop, wait=False)
		return protPop

	def test(self):
		protMHCII = self._runMHCIIPrediction()
		self._waitOutput(protMHCII, 'outputROIs', sleepTime=10)
		protPop = self._runMHCCoverage(protMHCII)
		self._waitOutput(protPop, 'outputROIs', sleepTime=10)
		assertHandle(self.assertIsNotNone, getattr(protPop, 'outputROIs', None))

class TestElliProPrediction(BaseTest):
	@classmethod
	def setUpClass(cls):
		super().setUpClass()
		cls.ds = DataSet.getDataSet('model_building_tutorial')
		setupTestProject(cls)

		cls._runImportPDB()
		cls._waitOutput(cls.protImportPDB, 'outputPdb', sleepTime=5)

	@classmethod
	def _runImportPDB(cls):
		cls.protImportPDB = cls.newProtocol(ProtImportPdb,
			inputPdbData=1, pdbFile=cls.ds.getFile('PDBx_mmCIF/5ni1.pdb'))
		cls.proj.launchProtocol(cls.protImportPDB, wait=False)

	def _runElliProPrediction(self):
		protElliPro = self.newProtocol(ProtElliProPrediction,
																	 rchains=True, chain_name='{"model": 0, "chain": "A", "residues": 92}')

		protElliPro.inputAtomStruct.set(self.protImportPDB)
		protElliPro.inputAtomStruct.setExtended('outputPdb')

		self.proj.launchProtocol(protElliPro, wait=False)
		return protElliPro

	def test(self):
		protElliPro = self._runElliProPrediction()
		self._waitOutput(protElliPro, 'outputROIs', sleepTime=10)
		assertHandle(self.assertIsNotNone, getattr(protElliPro, 'outputStructROIs', None))
		assertHandle(self.assertIsNotNone, getattr(protElliPro, 'outSequenceROIs_A', None))