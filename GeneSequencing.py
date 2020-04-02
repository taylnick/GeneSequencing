#!/usr/bin/python3

from which_pyqt import PYQT_VER

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
                        algo_output = self.banded_algorithm(sequences[i], sequences[j])
                    else:
                        algo_output = self.unrestricted_algorithm(sequences[i], sequences[j])

                    score, alignment1, alignment2 = algo_output

                    # Thank you, Tyler!
                    # print('\n' + str(i + 1) + ' vs ' + str(j + 1))
                    # print(alignment1[:100])
                    # print(alignment2[:100])
                    # print('Score: ' + str(score))

                    ###################################################################################################
                    s = {'align_cost': score, 'seqi_first100': alignment1, 'seqj_first100': alignment2}
                    table.item(i, j).setText('{}'.format(int(score) if score != math.inf else score))
                    table.update()
                jresults.append(s)
            results.append(jresults)
        return results

# Unrestricted algorithm: takes two strings to be compared and outputs the score, and corresponding alignment strings.
    # Author: Nathan Taylor CS 4412
    def unrestricted_algorithm(self, vert, hori):
        # Cost Table
        n = len(vert) + 1 if len(vert) < self.MaxCharactersToAlign else self.MaxCharactersToAlign + 1
        m = len(hori) + 1 if len(hori) < self.MaxCharactersToAlign else self.MaxCharactersToAlign + 1
        CT = [[0] * m for x in range(n)]
        # Back Pointer array
        BP = [['n'] * m for x in range(n)]
        vert_align = ""
        hori_align = ""

        # initialize values for row and column 0
        for each in range(n):
            CT[each][0] = each * INDEL
            BP[each][0] = 'a'
        for each in range(m):
            CT[0][each] = each * INDEL
            BP[0][each] = 'l'
        # Make sure the 0,0 spot is 'n'. The above loops reset it.
        BP[0][0] = 'n'
        # Fill in the rest of the table.
        for i in range(1, n):
            for j in range(1, m):
                match_or_sub = MATCH if vert[i - 1] == hori[j - 1] else SUB
                diag = match_or_sub + CT[i - 1][j - 1]
                left = INDEL + CT[i][j - 1]
                above = INDEL + CT[i - 1][j]
                min_num = min(left, above, diag)
                # Check what the minimum number is and assign the value to the cell along with an int
                # corresponding to the direction of the backpointer. d is diagonal, l is left, and a is above and n is none.

                if min_num == left:
                    CT[i][j] = left
                    BP[i][j] = 'l'
                elif min_num == above:
                    CT[i][j] = above
                    BP[i][j] = 'a'
                elif min_num == diag:
                    CT[i][j] = diag
                    BP[i][j] = 'd'

        # Create alignments from backpointer arrays.
        i = n - 1
        j = m - 1

        while i != 0 and j != 0:
            curr = BP[i][j]
            if curr == "d":
                i -= 1
                j -= 1
                vert_align = vert[i] + vert_align
                hori_align = hori[j] + hori_align
            elif curr == 'l':
                j -= 1
                vert_align = '-' + vert_align
                hori_align = hori[j] + hori_align
            elif curr == 'a':
                i -= 1
                hori_align = '-' + hori_align
                vert_align = vert[i] + vert_align
            elif curr == 'n':
                # unreachable, no solution.
                return [math.inf, 'No Alignment Possible', 'No Alignment Possible']

        # #return the cost to be populated in the GUI.
        # vert_align and hori_align are needed to show what the end result is.
        cost = CT[-1][-1]
        return [cost, vert_align, hori_align]

# Banded n x k algorithm. Takes two strings and outputs the score with corresponding alignment strings.
    # Author: Nathan Taylor CS 4412
    def banded_algorithm(self, vert, hori):
        # initialize n and m to the proper lengths.
        n = len(vert) + 1 if len(vert) < self.MaxCharactersToAlign else self.MaxCharactersToAlign + 1
        m = len(hori) + 1 if len(hori) < self.MaxCharactersToAlign else self.MaxCharactersToAlign + 1
        k = 2 * MAXINDELS + 1

        # Cost Table
        CT = [[math.inf] * k for x in range(n)]
        # Back Pointer array
        BP = [['n'] * k for x in range(n)]
        vert_align = ""
        hori_align = ""

        # initialize values for row 0 and column 0 down to maxindel.
        # Using Maxindel accounts for k
        for each in range(MAXINDELS + 1):
            CT[each][0] = each * INDEL
            BP[each][0] = 'a'
        for each in range(MAXINDELS + 1):
            CT[0][each] = each * INDEL
            BP[0][each] = 'l'
        # Make sure the 0,0 spot is 'n'. The above loops reset it.
        BP[0][0] = 'n'

        # Fill in the rest of the table.

        # Top part.
        for i in range(1, MAXINDELS + 1):
            for j in range(1, MAXINDELS + i + 1):
                match_or_sub = MATCH if vert[i - 1] == hori[j - 1] else SUB
                diag = match_or_sub + CT[i - 1][j - 1]
                left = INDEL + CT[i][j - 1]
                above = INDEL + CT[i - 1][j]
                min_num = min(left, above, diag)
                # Check what the minimum number is and assign the value to the cell along with an int
                # corresponding to the direction of the backpointer. d is diagonal, l is left, and a is above and n is none.

                if min_num == left:
                    CT[i][j] = left
                    BP[i][j] = 'l'
                elif min_num == above:
                    CT[i][j] = above
                    BP[i][j] = 'a'
                elif min_num == diag:
                    CT[i][j] = diag
                    BP[i][j] = 'd'

        # Middle part:
        for i in range(MAXINDELS + 1, n - MAXINDELS):
            offset = i - MAXINDELS - 1
            for j in range(k):
                if j == 0:
                    left = math.inf
                    above = INDEL + CT[i - 1][j + 1]
                elif j < k - 1:
                    left = INDEL + CT[i][j - 1]
                    above = INDEL + CT[i - 1][j + 1]
                else:
                    left = INDEL + CT[i][j - 1]
                    above = math.inf
                match_or_sub = MATCH if vert[i - 1] == hori[j + offset] else SUB
                diag = match_or_sub + CT[i - 1][j]
                min_num = min(left, above, diag)

                if min_num == left:
                    CT[i][j] = left
                    BP[i][j] = 'l'
                elif min_num == above:
                    CT[i][j] = above
                    BP[i][j] = 'a'
                elif min_num == diag:
                    CT[i][j] = diag
                    BP[i][j] = 'd'

        # Bottom Part
        for i in range(n - MAXINDELS, n):
            #tail = i - (n - MAXINDELS)
            offset = i - MAXINDELS - 1
            for j in range(k):
                # n is length of the word along the vertical axis
                # m is length of word along horizontal axis
                if m - (j + offset) <= 1:
                    continue
                if j == 0:
                    left = math.inf
                    above = INDEL + CT[i - 1][j + 1]
                elif j < k - 1:
                    left = INDEL + CT[i][j - 1]
                    above = INDEL + CT[i - 1][j + 1]
                else:
                    left = INDEL + CT[i][j - 1]
                    above = math.inf
                match_or_sub = MATCH if vert[i - 1] == hori[j + offset] else SUB
                diag = match_or_sub + CT[i - 1][j]
                min_num = min(left, above, diag)
                # Check what the minimum number is and assign the value to the cell along with an int corresponding
                # to the direction of the backpointer. d is diagonal, l is left, and a is above and n is none.

                if min_num == left:
                    CT[i][j] = left
                    BP[i][j] = 'l'
                elif min_num == above:
                    CT[i][j] = above
                    BP[i][j] = 'a'
                elif min_num == diag:
                    CT[i][j] = diag
                    BP[i][j] = 'd'

        # Create alignments from backpointer arrays.
        i = n - 1
        j = m - 1
        if abs(i - j) >= 4:
            return [math.inf, 'No Alignment Possible', 'No Alignment Possible']

        # Find the last non null entry to use as a beginning point for back propogation.
        x = n - 1
        y = k - 1
        direction = BP[x][y]
        while direction == 'n':
            y -= 1
            direction = BP[x][y]
        while x > MAXINDELS:
            direction = BP[x][y]
            if direction == "d":
                i -= 1
                j -= 1
                vert_align = vert[i] + vert_align
                hori_align = hori[j] + hori_align
                x -= 1

            elif direction == 'l':
                j -= 1
                vert_align = '-' + vert_align
                hori_align = hori[j] + hori_align
                y -= 1
            elif direction == 'a':
                i -= 1
                hori_align = '-' + hori_align
                vert_align = vert[i] + hori_align
                x -= 1
                y += 1
            elif direction == 'n':# unreachable, no solution.
                return [math.inf, 'No Alignment Possible', 'No Alignment Possible']
        i -= 1
        y -= 1
        while i != 0 and y != 0:
            curr = BP[i][y]
            if curr == "d":
                i -= 1
                y -= 1
                vert_align = vert[i] + vert_align
                hori_align = hori[y] + hori_align
            elif curr == 'l':
                y -= 1
                vert_align = '-' + vert_align
                hori_align = hori[y] + hori_align
            elif curr == 'a':
                i -= 1
                hori_align = '-' + hori_align
                vert_align = vert[i] + vert_align
            elif curr == 'n':
                # unreachable, no solution.
                return [math.inf, 'No Alignment Possible', 'No Alignment Possible']

        # return the cost to be populated in the GUI.
        # Look for the last non-infinite entry.
        last = -1
        cost = CT[-1][last]
        while cost == math.inf:
            last -= 1
            cost = CT[-1][last]

        return [cost, vert_align, hori_align]
