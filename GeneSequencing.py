#!/usr/bin/python3

from which_pyqt import PYQT_VER
from numpy import numarray

if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1


class GeneSequencing:

    def __init__(self):
        self.banded = bool
        self.MaxCharactersToAlign = int

    # This is the method called by the GUI.  _sequences_ is a list of the ten sequences, _table_ is a
    # handle to the GUI so it can be updated as you find results, _banded_ is a boolean that tells
    # you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
    # how many base pairs to use in computing the alignment

    def align(self, sequences, table, banded, align_length):
        self.banded = banded
        self.MaxCharactersToAlign = align_length
        results = []

        for i in range(len(sequences)):
            jresults = []
            for j in range(len(sequences)):

                if j < i:
                    s = {}
                else:
                    ###################################################################################################
                    # your code should replace these three statements and populate the three variables:
                    # score, alignment1 and alignment2
                    algo_output = []
                    if banded:
                        algo_output = self.banded_algorithm(sequences.pop(), sequences[j])
                    else:
                        algo_output = self.unrestricted_algorithm(sequences[i], sequences[j])

                    score = algo_output[0]
                    alignment1 = 'abc-easy  DEBUG:(seq{}, {} chars,align_len={}{})'.format(i + 1,
                                                                                           len(sequences[i]),
                                                                                           align_length,
                                                                                           ',BANDED' if banded else '')
                    alignment2 = 'as-123--  DEBUG:(seq{}, {} chars,align_len={}{})'.format(j + 1,
                                                                                           len(sequences[j]),
                                                                                           align_length,
                                                                                           ',BANDED' if banded else '')
                    ###################################################################################################
                    s = {'align_cost': score, 'seqi_first100': alignment1, 'seqj_first100': alignment2}
                    table.item(i, j).setText('{}'.format(int(score) if score != math.inf else score))
                    table.update()
                jresults.append(s)
            results.append(jresults)
        return results

    # TODO: Implement unrestricted_algorithm
    def unrestricted_algorithm(self, vert, hori):
        # Cost Table
        n = len(vert) if len(vert) < self.MaxCharactersToAlign else self.MaxCharactersToAlign
        m = len(hori) if len(hori) < self.MaxCharactersToAlign else self.MaxCharactersToAlign
        CT = [[0] * m for x in range(n)]

        vert_align = ""
        hori_align = ""

        # initialize values for row and column 0
        for each in range(n):
            CT[each][0] = each
        for each in range(m):
            CT[0][each] = each

        # Fill in the rest of the table.
        for i in range(1, n):
            for j in range(1, m):
                leftDiag = self.match_sub(vert[i], hori[j]) + CT[i - 1][j - 1]
                CT[i][j] = min(leftDiag, INDEL + CT[i][j - 1],
                               INDEL + CT[i - 1][j])

        #return the cost to be populated in the GUI.
        # vert_align and hori_align are needed to show what the end result is.
        cost = CT[-1][-1]
        return [cost, vert_align, hori_align]

    # function to check if two characters are identical or if a substitution is necessary.
    # returns an int.
    def match_sub(self, char1, char2):
        return MATCH if char1 == char2 else SUB

    # TODO: Implement banded_algorithm
    def banded_algorithm(self, vert, hori):
        cost = int
        vert_align = str
        hori_align = str

        return [cost, vert_align, hori_align]
