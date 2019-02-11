#!/usr/bin/env python

import numpy as np

#penalty box
match    = 3
mismatch_transit = -1.5
mismatch_transver = -3
gap      = -2

#sequences. seq2 should always be larger than seq1
seq1 = 'TGTTACGG'
seq2 = 'GGTTGACTA'

#scoring matrix size
rows = len(seq1) + 1
cols = len(seq2) + 1

def create_score_matrix(rows, cols):
    '''
    Create a matrix of scores representing trial alignments of the two sequences.
    Sequence alignment can be treated as a graph search problem. This function
    creates a graph (2D matrix) of scores, which are based on trial alignments
    of different base pairs. The path with the highest cummulative score is the
    best alignment.
    '''
    score_matrix = [[0 for col in range(cols)] for row in range(rows)]
    # Fill the scoring matrix.
    max_score = 0
    max_pos   = None    # The row and columbn of the highest score in matrix.
    for i in range(1, rows):
        for j in range(1, cols):
            score = calc_score(score_matrix, i, j)
            if score > max_score:
                max_score = score
                max_pos   = (i, j)
            score_matrix[i][j] = score
    assert max_pos is not None, 'the x, y position with the highest score was not found'
    return score_matrix, max_pos

def calc_score(matrix, x, y):
    #Calculate score for a given x, y position in the scoring matrix.
    #The score is based on the up, left, and upper-left diagnol neighbors
    if seq1[x - 1] == seq2[y - 1]:
        similarity = match
    elif seq1[x - 1] == 'A' and seq2[y - 1] == 'G':
        similarity = mismatch_transit
    elif seq1[x - 1] == 'C' and seq2[y - 1] == 'T':
        similarity = mismatch_transit
    else:
        similarity = mismatch_transver
    diag_score = matrix[x - 1][y - 1] + similarity
    up_score   = matrix[x - 1][y] + gap
    left_score = matrix[x][y - 1] + gap
    return max(0, diag_score, up_score, left_score)

def print_matrix(matrix):
    '''
    Print the scoring matrix.
    ex:
    0   0   0   0   0   0
    0   2   1   2   1   2
    0   1   1   1   1   1
    0   0   3   2   3   2
    0   2   2   5   4   5
    0   1   4   4   7   6
    '''
    print(np.matrix(matrix).T)

#add your function(s) to find a solution here.

score_matrix, start_pos = create_score_matrix(rows, cols)

def find_path(matrix, start_pos):
    '''
    print a string of a path with arrows (i1j1 -> i2j2 -> i3j3). From the biggest number to 0, which is the starting point.
    If the three number around the biggest number are same, you can choose a random one to go.
    '''
    x, y = start_pos[0], start_pos[1]
    curr_point = matrix[x][y]
    loc = []
    loc.append((x,y))
    while curr_point != 0:
        potential = []
        x_up, y_up = x, y-1
        x_left, y_left = x-1, y
        x_dia, y_dia = x-1, y-1
        potential.append((matrix[x_up][y_up], x_up, y_up))
        potential.append((matrix[x_left][y_left], x_left, y_left))
        potential.append((matrix[x_dia][y_dia], x_dia, y_dia))
        target = max(potential, key = lambda x:x[0])
        curr_point = target[0]
        x, y = target[1], target[2]
        loc.append((x,y))
    result = 'The path is: '
    for coord in loc:
        if coord == loc[-1]:
            result += str(coord)
        else:
            result += str(coord) + ' -> '
    print(result)

#end of your function(s)

if __name__ == '__main__':
    #my main
    score_matrix, start_pos = create_score_matrix(rows, cols)
    print_matrix(score_matrix)
    # find_path(score_matrix, start_pos)

