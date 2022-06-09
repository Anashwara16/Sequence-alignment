# Sequence-alignment

Implemented Needleman–Wunsch algorithm and Hirschberg's algorithm for Optimal Sequence Alignment Problem in Python.

**Analysis & Observations:**

**Observations: **
    1. There can be multiple alignments for a given set of input strings (with the minimum cost).
    2. Basic version of the algorithm (using Dynamic Programming) generates a minimal cost alignment in lesser time when compared to Memory efficient version but uses more space.
    3. Memory efficient version of the algorithm (using Divide and Conquer) is expected to generate the sequence alignment with lesser memory utilization when compared to Dynamic programming version. However, both the algorithms take similar amount of space when the input size is small (<400 in our plots).

**Analysis: **
    1. The basic version takes O(mn) memory and has a running time of O(mn), because at worst, a constant time is spent on each element of the O(mn) array. 
    2. The memory efficient version takes a constant time factor more than the basic version and using the masters theorem, we can compute the time complexity to be O(mn) for this Divide and Conquer algorithm. 
    3. This memory efficient algorithm takes O(m+n) space.
