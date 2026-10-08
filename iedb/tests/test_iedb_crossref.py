# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Enzo Sierra (enzogael57@gmail.com)
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
import unittest

from pyworkflow.tests import setupTestProject, BaseTest

from pwem.protocols import ProtImportSequence
from pwchem.protocols import ProtDefineSeqROI

from .. import Plugin as iedbPlugin
from ..constants import IEDB_BCELL_REFERENCE_PATH
from ..protocols import ProtIEDBCrossref


@unittest.skipUnless(
    os.path.isfile(iedbPlugin.getVar(IEDB_BCELL_REFERENCE_PATH) or ''),
    'IEDB_BCELL_REFERENCE_PATH not configured (manually-downloaded/filtered IEDB B-cell reference)')
class TestIEDBCrossref(BaseTest):
    # 'TMSPISIADRDFGIDIPNIPQ' -- a real documented neutralizing linear
    # epitope (Orientia tsutsugamushi, PMID 10924780), confirmed present
    # verbatim in the configured IEDB_BCELL_REFERENCE_PATH reference CSV.
    #
    # NO_MATCH_SEGMENT is deliberately NOT a homopolymer (e.g. 'AAAA...'):
    # a real run against the actual reference file found a genuine
    # 2-residue 'AA' epitope entry, whose 2 < minOverlap(=6) exemption
    # (reference epitopes SHORTER than minOverlap only need their own full
    # length matched, see queryIedbCrossref's docstring) makes ANY run of
    # >=2 identical residues a spurious "match" -- confirmed by reading the
    # real 'extra/iedb_crossref.csv' from that run, not assumed. The 20
    # standard amino acids used once each, in this fixed order, has no
    # 6-residue substring anywhere in the real reference file (grep-
    # verified) -- a genuine negative control.
    NAME = 'IEDB_CROSSREF_TEST_SEQ'
    MATCH_SEGMENT = 'TMSPISIADRDFGIDIPNIPQ'
    NO_MATCH_SEGMENT = 'ACDEFGHIKLMNPQRSTVWY'
    AMINOACIDSSEQ = 'GGGGG' + MATCH_SEGMENT + 'GGGGG' + NO_MATCH_SEGMENT

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        setupTestProject(cls)

        cls._runImportSeq()
        cls._waitOutput(cls.protImportSeq, 'outputSequence', sleepTime=5)

    @classmethod
    def _runImportSeq(cls):
        kwargs = {
            'inputSequenceName': cls.NAME,
            'inputSequenceDescription': cls.NAME,
            'inputRawSequence': cls.AMINOACIDSSEQ,
        }
        cls.protImportSeq = cls.newProtocol(ProtImportSequence, **kwargs)
        cls.proj.launchProtocol(cls.protImportSeq, wait=False)

    @classmethod
    def _runDefSeqROIs(cls, inProt):
        matchStart = cls.AMINOACIDSSEQ.index(cls.MATCH_SEGMENT) + 1
        matchEnd = matchStart + len(cls.MATCH_SEGMENT) - 1
        noMatchStart, noMatchEnd = matchEnd + 6, len(cls.AMINOACIDSSEQ)
        inROIs = (
            '1) Residues: {{"index": "{}-{}", "residues": "{}", "desc": "None"}}\n'
            '2) Residues: {{"index": "{}-{}", "residues": "{}", "desc": "None"}}'
        ).format(
            matchStart, matchEnd, cls.MATCH_SEGMENT,
            noMatchStart, noMatchEnd, cls.AMINOACIDSSEQ[noMatchStart - 1:noMatchEnd],
        )
        protDefSeqROIs = cls.newProtocol(ProtDefineSeqROI, chooseInput=0, inROIs=inROIs)
        protDefSeqROIs.inputSequence.set(inProt)
        protDefSeqROIs.inputSequence.setExtended('outputSequence')
        cls.proj.launchProtocol(protDefSeqROIs, wait=False)
        return protDefSeqROIs

    def _runIEDBCrossref(self, protROIs):
        protCrossref = self.newProtocol(ProtIEDBCrossref)
        protCrossref.inputROIs.set(protROIs)
        protCrossref.inputROIs.setExtended('outputROIs')
        self.proj.launchProtocol(protCrossref, wait=False)
        return protCrossref

    def test(self):
        protROIs = self._runDefSeqROIs(inProt=self.protImportSeq)
        self._waitOutput(protROIs, 'outputROIs', sleepTime=5)

        protCrossref = self._runIEDBCrossref(protROIs)
        self._waitOutput(protCrossref, 'outputROIs', sleepTime=5)

        outROIs = getattr(protCrossref, 'outputROIs', None)
        self.assertIsNotNone(outROIs)

        bestLenBySeq = {roi.getROISequence(): roi._iedbBestMatchLength.get() for roi in outROIs}
        bestOrgBySeq = {roi.getROISequence(): roi._iedbBestOrganism.get() for roi in outROIs}
        self.assertEqual(len(bestLenBySeq), 2)

        # The real documented epitope: full-length match, correct source
        # organism.
        self.assertEqual(bestLenBySeq[self.MATCH_SEGMENT], len(self.MATCH_SEGMENT))
        self.assertEqual(bestOrgBySeq[self.MATCH_SEGMENT], 'Orientia tsutsugamushi str. Boryong')

        # NO_MATCH_SEGMENT: no 6+ residue substring of it exists anywhere in
        # the real reference file (grep-verified, see class docstring), but
        # the reference file does contain some genuine epitope entries
        # SHORTER than minOverlap (e.g. a real 2-residue 'ST' entry) --
        # queryIedbCrossref's documented behavior is that those still count
        # as a match on their own full (short) length. So the correct
        # invariant here is "no match anywhere near the real segment's
        # length", not "zero matches".
        noMatchSeq = next(seq for seq in bestLenBySeq if seq != self.MATCH_SEGMENT)
        self.assertLess(bestLenBySeq[noMatchSeq], 6)
