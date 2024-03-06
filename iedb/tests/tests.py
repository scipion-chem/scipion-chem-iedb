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

from pwem.protocols import ProtImportSequence

from pwchem.utils import assertHandle

from ..protocols import ProtBepiPredPrediction, ProtMHCIPrediction


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