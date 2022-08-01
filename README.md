# Sequence-alignment

Implemented Needleman–Wunsch algorithm and Hirschberg's algorithm for Optimal Sequence Alignment Problem in Python.

Analysis & Observations:

Analysis: 

The basic version takes O(mn) memory and has a running time of O(mn), because at worst, a constant time is spent on each element of the O(mn) array. The memory efficient version takes a constant time factor more than the basic version. Using the masters theorem, we can compute the time complexity to be O(mn) for this Divide and Conquer algorithm. This memory efficient algorithm takes O(m+n) space.
    
Observations:

There can be multiple alignments for a given set of input strings (with the minimum cost). Basic version of the algorithm (using Dynamic Programming) generates a minimal cost alignment in lesser time when compared to memory efficient version but uses more space. Memory efficient version of the algorithm (using Divide and Conquer) is expected to generate the sequence alignment with lesser memory utilization when compared to Dynamic programming version. However, both the algorithms take similar amount of space when the input size is small (<400 in plots).
 
![Alt text](https://github.com/Anashwara16/Sequence-alignment/blob/main/CPUPlot.png?raw=true)


![Alt text](https://github.com/Anashwara16/Sequence-alignment/blob/main/MemoryPlot.png?raw=true)

