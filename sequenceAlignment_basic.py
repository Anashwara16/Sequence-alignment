#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import pdb
import time
import os
import psutil
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



def profileDP(X, Y, report=False):
    process = psutil.Process(os.getpid())
    startTime = time.process_time()
    tracemalloc.start()
    result = SequenceAlignmentDP(X, Y)
    executionTime = (time.process_time() - startTime)
    current,peak=tracemalloc.get_traced_memory()
    tracemalloc.stop()
    memoryUsed = float(peak/5)
    output(result, executionTime, memoryUsed, report)
    # return (len(X) + len(Y), executionTime, memoryUsed)


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
        print(float(memoryUsed), file=of)
        print('\n')


if __name__ == '__main__':

    X, Y = stringGenerator(sys.argv[1])
    profileDP(X, Y)
