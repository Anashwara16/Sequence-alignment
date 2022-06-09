#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import enum
import pdb
import time
import os
import psutil
# from util import output
import tracemalloc

gap_penalty = 30

mismatch_cost = {'AA': 0, 'AC': 110, 'AG': 48, 'AT': 94,
                 'CA': 110, 'CC': 0, 'CG': 118, 'CT': 48,
                 'GA': 48, 'GC': 118, 'GG': 0, 'GT': 110,
                 'TA': 94, 'TC': 48, 'TG': 110, 'TT': 0}

def stringGenerator(inputPath):
    lines = open(inputPath, 'r').read().splitlines()
    for line in lines:
        try:
            index = int(line)
            if index >= 0 and index < len(generatedString):
                generatedString = generatedString[0:index + 1] + \
                generatedString + generatedString[index + 1:]

        except:
            if lines.index(line) != 0:
                firstString = generatedString
            baseString = line.upper()
            generatedString = baseString
    # print("##### In string generator####")
    # print(firstString, generatedString, sep = '\n')
    # print("##### Exiting string generator####")
    return firstString, generatedString


def SequenceAlignmentDP(X, Y):
    m = len(X)
    n = len(Y)
    values = [[0 for x in range(n+1)] for y in range(m+1)]
    for i in range(m+1):
        values[i][0] = i * gap_penalty
    for j in range(n+1):
        values[0][j] = j * gap_penalty
    for j in range(1, n+1):
        for i in range(1, m+1):
            values[i][j] = min((values[i-1][j-1] + mismatch_cost[X[i-1] + Y[j-1]]),
                               (values[i-1][j] + gap_penalty),
                               (values[i][j-1] + gap_penalty))

    # print(values)
    i = m
    j = n
    align1 = ''
    align2 = ''

    while i > 0 and j > 0:
        if values[i][j] == values[i-1][j-1] + mismatch_cost[X[i-1] + Y[j-1]]:
            align1 = X[i-1] + align1
            align2 = Y[j-1] + align2
            i -= 1
            j -= 1
        elif values[i][j] == values[i-1][j] + gap_penalty:
            align1 = X[i-1] + align1
            align2 = '_' + align2
            i -= 1
        elif values[i][j] == values[i][j-1] + gap_penalty:
            align1 = '_' + align1
            align2 = Y[j-1] + align2
            j -= 1

    if i > 0:
        align1 = X[0:i] + align1
        align2 = '_'*i + align2
    if j > 0:
        align2 = Y[0:j] + align2
        align1 = '_'*j + align1
    # print("######## INside DP func")
    # print(align1)
    # print(align2)
    # print(Test_Cost(align1, align2))
    # print("######## Exiting DP func")
    return[align1, align2,values[m][n]]


class AlignType(enum.Enum):
    ForwardAlignment = 1
    BackwardAlignment = 2


def spaceEfficientAlign(X, Y, align):
    m, n = len(X), len(Y)

    # Array B[0...m, 0...1]
    B = [[0] * 2 for i in range(m + 1)]

    # Initialize B[i,0]=iδ for each i (just as in column 0 of A)
    for i in range(m+1):
        B[i][0] = gap_penalty * i

    for j in range(1, n + 1):
        # B[0,1]=jδ (since this corresponds to entry A[0,j])
        B[0][1] = j * gap_penalty
        for i in range(1, m + 1):
            if (align == AlignType.ForwardAlignment):
                mismatchCost = mismatch_cost[X[i - 1] + Y[j - 1]]
            if (align == AlignType.BackwardAlignment):
                mismatchCost = mismatch_cost[X[m - i] + Y[n - j]]
            B[i][1] = min(B[i-1][0] + mismatchCost,
                          B[i-1][1] + gap_penalty,
                          B[i][0] + gap_penalty)

        # Move column 1 of B to column 0 to make room for next iteration:
        for i in range(m+1):
            B[i][0] = B[i][1]

    # Return 1st column of B
    return [row[1] for row in B]


def SequenceAlignmentDC(X, Y):
    # Let m be the number of symbols in X
    m = len(X)

    # Let n be the number of symbols in Y
    n = len(Y)

    # Base Case:
    if n <= 2 or m <= 2:
        # For n,m < 2 use the dynamic programming sequence alignment algorithm
        return SequenceAlignmentDP(X, Y)
    midY = int(n / 2)

    # Call Space-Efficient-Alignment( X, Y[ 1 : n/2 ] )
    fAlign = spaceEfficientAlign(X, Y[:midY], AlignType.ForwardAlignment)

    # Call Backward-Space-Efficient-Alignment(X, Y[ n/2 + 1 : n ])
    bAlign = spaceEfficientAlign(X, Y[midY:], AlignType.BackwardAlignment)

    # Let q be the index minimizing f(q, n/2) + g(q, n/2)
    alignSum = [fAlign[j] + bAlign[m-j] for j in range(m + 1)]
    q = alignSum.index(min(alignSum))

    # Clear temporary lists.
    fAlign, bAlign, alignSum = [], [], []

    # Divide-and-Conquer-Alignment(X[ 1 : q ], Y[ 1 : n/2 ])
    left = SequenceAlignmentDC(X[:q], Y[:midY])
    # Divide-and-Conquer-Alignment(X[ q + 1 : n ],Y[ n/2 + 1 : n ])
    right = SequenceAlignmentDC(X[q:], Y[midY:])

    # Return result in format: [1st alignment, 2nd alignment, similarity]
    return [left[r] + right[r] for r in range(3)]

def output(result, executionTime, memoryUsed, report=False):
    # align1 = result[0].lstrip('_')
    # align2 = result[1].lstrip('_')
    align1=result[0]
    align2=result[1]
    # if report is False:
    #     print(align1)
    #     print(align2)
    #     print(result[2])
    with open('output.txt', 'w') as of:
        print(align1[:50] + " " + align1[-50:], file=of)
        print(align2[:50] + " " + align2[-50:], file=of)
        print(result[2], file=of)
        print(executionTime, file=of)
        print(memoryUsed, file=of)
        print('\n')

def profileDC(X, Y, report=False):
    process = psutil.Process(os.getpid())
    startTime = time.process_time()
    tracemalloc.start()
    result = SequenceAlignmentDC(X, Y)
    executionTime = (time.process_time() - startTime)
    current,peak=tracemalloc.get_traced_memory()
    tracemalloc.stop()
    memoryUsed = float(peak)
    output(result, executionTime, memoryUsed)
    # return (len(X) + len(Y), executionTime, memoryUsed)




if __name__ == '__main__':
    X, Y = stringGenerator(sys.argv[1])
    profileDC(X, Y)
