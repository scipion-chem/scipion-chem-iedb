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

"""
This protocol is used to cross-reference peptide candidates against known
B-cell epitopes with a documented protective or neutralizing effect, for
any pathogen, using a locally-filtered subset of the IEDB bulk export.

Generalizes the HIV-specific LANL/CATNAP cross-reference to any studied
pathogen: same matching mechanism (longest-common-substring overlap, no
structural alignment or organism pre-filtering -- the sequence match
itself acts as the implicit organism filter). It is self-contained: the
matching helper is defined below rather than imported from another plugin,
so this protocol carries no cross-plugin dependency.

IEDB's bulk export terms of use do not clearly permit redistributing even
a filtered subset, so -- same treatment as LANL/CATNAP -- this stays a
locally-downloaded reference file, never auto-fetched.
"""

import os
from typing import List

import pandas as pd
from pwchem.objects import SetOfSequenceROIs
from pwem.protocols import EMProtocol
from pyworkflow.object import Integer, String
from pyworkflow.protocol import params
from pyworkflow.utils import Message

from .. import Plugin as iedbPlugin
from ..constants import IEDB_BCELL_REFERENCE_PATH

_OUTPUT_COLUMNS = [
    'sequence', 'epitope_sequence', 'match_length', 'source_organism',
    'response_measured', 'qualitative_measure', 'method', 'host', 'pmid',
]

_RAW_CSV_COLUMNS = [
    'epitope_sequence', 'source_organism', 'response_measured',
    'qualitative_measure', 'method', 'host', 'pmid',
]

DEFAULT_MIN_OVERLAP = 6


def longestCommonSubstringLen(a: str, b: str) -> int:
    """Length of the longest common substring between 'a' and 'b' (DP O(len(a)*len(b)))."""
    if not a or not b:
        return 0
    prev = [0] * (len(b) + 1)
    best = 0
    for i in range(1, len(a) + 1):
        curr = [0] * (len(b) + 1)
        for j in range(1, len(b) + 1):
            if a[i - 1] == b[j - 1]:
                curr[j] = prev[j - 1] + 1
                best = max(best, curr[j])
        prev = curr
    return best


class IEDBParseError(Exception):
    """The IEDB reference file does not match the expected format."""


def loadIedbEpitopes(iedbCsvPath: str) -> pd.DataFrame:
    """Load the already-filtered IEDB B-cell protective/neutralizing epitope subset."""
    df = pd.read_csv(iedbCsvPath, dtype=str)
    missing = [c for c in _RAW_CSV_COLUMNS if c not in df.columns]
    if missing:
        raise IEDBParseError(f"'{iedbCsvPath}' does not have the expected columns: missing {missing}.")
    df['epitope_sequence'] = df['epitope_sequence'].str.strip().str.upper()
    return df


def queryIedbCrossref(
    sequences: List[str], iedbCsvPath: str, minOverlap: int = DEFAULT_MIN_OVERLAP,
) -> pd.DataFrame:
    """Cross-reference 'sequences' against IEDB epitopes with a documented protective/neutralizing effect.

    Does NOT pre-filter by organism: the sequence match itself acts as the
    implicit filter (same principle as ProtLANLCATNAPCrossref).

    Args:
        sequences: Candidate peptides/sequences to evaluate.
        iedbCsvPath: Path to the locally-filtered IEDB B-cell reference CSV.
        minOverlap: Minimum substring overlap length to report a match.
            For reference epitopes SHORTER than this, the full epitope
            match is required instead (never a looser threshold than the
            epitope itself).

    Returns:
        DataFrame with one row per (candidate, reference epitope) pair
        that overlaps enough, columns _OUTPUT_COLUMNS. Empty if no match
        or 'sequences' is empty.
    """
    if not sequences:
        return pd.DataFrame(columns=_OUTPUT_COLUMNS)

    iedbDf = loadIedbEpitopes(iedbCsvPath)
    if iedbDf.empty:
        return pd.DataFrame(columns=_OUTPUT_COLUMNS)

    rows = []
    for seq in sequences:
        seqUpper = seq.upper()
        for ref in iedbDf.itertuples(index=False):
            requiredOverlap = min(minOverlap, len(ref.epitope_sequence))
            matchLen = longestCommonSubstringLen(seqUpper, ref.epitope_sequence)
            if matchLen < requiredOverlap:
                continue

            rows.append({
                'sequence': seq,
                'epitope_sequence': ref.epitope_sequence,
                'match_length': matchLen,
                'source_organism': ref.source_organism,
                'response_measured': ref.response_measured,
                'qualitative_measure': ref.qualitative_measure,
                'method': ref.method,
                'host': ref.host,
                'pmid': ref.pmid,
            })

    return pd.DataFrame(rows, columns=_OUTPUT_COLUMNS) if rows else pd.DataFrame(columns=_OUTPUT_COLUMNS)


class ProtIEDBCrossref(EMProtocol):
    """
    AI Generated:

    Cross-references every input ROI's peptide against IEDB B-cell
    epitopes with a documented protective/neutralizing effect (any
    pathogen, not HIV-specific -- see pwchem core's ProtLANLCATNAPCrossref
    for that), and annotates (does NOT filter) each ROI with a summary of
    its best match. Purely informative: harmless (zero matches) when the
    input pathogen has no documented epitope in the local reference subset.

    Output
    ------
    outputROIs: the same SetOfSequenceROIs as the input, annotated with
    '_iedbMatchCount' (int, total matching reference epitopes),
    '_iedbBestMatchLength' (int, the longest match), '_iedbBestOrganism'
    (str, source organism of that best match). The full detail (every
    match, not just the best one) is persisted to
    'extra/iedb_crossref.csv'.
    """

    _label = 'iedb b-cell crossref'

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputROIs', params.PointerParam, pointerClass='SetOfSequenceROIs',
                       label='Sequence ROIs: ',
                       help='Peptide candidates to cross-reference against IEDB documented '
                            'protective/neutralizing epitopes.')
        form.addParam('minOverlap', params.IntParam, default=DEFAULT_MIN_OVERLAP,
                       label='Min. substring overlap (aa): ',
                       help='Minimum substring overlap to report a match. Reference epitopes '
                            'shorter than this still require their own full length as the match.')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.crossrefStep)
        self._insertFunctionStep(self.createOutputStep)

    # ---------------------------------- Steps -----------------------------------

    def _getCrossrefPath(self):
        return self._getExtraPath('iedb_crossref.csv')

    def _getRois(self):
        # Iterating a Scipion SetOfXXX reuses the same Python object per row
        # (the underlying sqlite cursor): each item must be cloned when
        # materialized into a list, or all N references end up pointing to
        # the cursor's last state.
        return [roi.clone() for roi in self.inputROIs.get()]

    def crossrefStep(self):
        rois = self._getRois()
        sequences = [roi.getROISequence() for roi in rois]
        if not sequences:
            return

        crossrefDf = queryIedbCrossref(
            sequences, iedbCsvPath=iedbPlugin.getVar(IEDB_BCELL_REFERENCE_PATH),
            minOverlap=self.minOverlap.get(),
        )
        crossrefDf.to_csv(self._getCrossrefPath(), index=False)

    def createOutputStep(self):
        rois = self._getRois()
        crossrefDf = pd.read_csv(self._getCrossrefPath()) if os.path.isfile(self._getCrossrefPath()) else pd.DataFrame()

        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
        for roi in rois:
            matches = crossrefDf[crossrefDf['sequence'] == roi.getROISequence()] if not crossrefDf.empty else crossrefDf

            roi._iedbMatchCount = Integer(len(matches))
            if len(matches):
                best = matches.sort_values('match_length', ascending=False).iloc[0]
                roi._iedbBestMatchLength = Integer(int(best['match_length']))
                roi._iedbBestOrganism = String(best['source_organism'])
            else:
                roi._iedbBestMatchLength = Integer(0)
                roi._iedbBestOrganism = String('')
            outROIs.append(roi)

        if len(outROIs) > 0:
            self._defineOutputs(outputROIs=outROIs)
            self._defineSourceRelation(self.inputROIs, outROIs)

    # ---------------------------------- Validation -------------------------------

    def _validate(self):
        errors = []
        iedbPath = iedbPlugin.getVar(IEDB_BCELL_REFERENCE_PATH)
        if not iedbPath or not os.path.isfile(iedbPath):
            errors.append(
                f"IEDB_BCELL_REFERENCE_PATH is not set or does not exist: '{iedbPath}'. Download the "
                "IEDB B-cell bulk export from https://www.iedb.org/database_export_v3.php, filter to "
                "linear peptides with a documented protective/neutralizing effect, and set "
                "IEDB_BCELL_REFERENCE_PATH in scipion.conf."
            )
        return errors

    def _summary(self):
        summary = []
        if self.isFinished():
            outROIs = getattr(self, 'outputROIs', None)
            if outROIs is not None:
                nMatch = sum(1 for roi in outROIs if roi._iedbMatchCount.get() > 0)
                summary.append(f'{nMatch}/{len(outROIs)} candidate(s) match >= 1 documented '
                               f'protective/neutralizing IEDB epitope.')
        return summary
