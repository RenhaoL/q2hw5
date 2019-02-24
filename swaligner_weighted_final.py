#!/usr/bin/env python

import numpy as np

#penalty box
match    = 4
mismatch_transit = -2
mismatch_transver = -4
gap      = -2

#sequences. seq2 should always be larger than seq1

# seq1 = 'GACTTAC'
# seq2 = 'CGTGAATTCAT'
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
    max_pos   = None    # The row and column of the highest score in matrix.
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
    elif seq1[x - 1] == 'A' and seq2[y - 1] == 'G' or seq1[x - 1] == 'G' and seq2[y - 1] == 'A': # give a different score to transition
        similarity = mismatch_transit
    elif seq1[x - 1] == 'C' and seq2[y - 1] == 'T' or seq1[x - 1] == 'T' and seq2[y - 1] == 'C': # give a different score to transition
        similarity = mismatch_transit
    else: # give a different score to transversion
        similarity = mismatch_transver
    diag_score = matrix[x - 1][y - 1] + similarity
    up_score   = matrix[x - 1][y] + gap
    left_score = matrix[x][y - 1] + gap
    return max(diag_score, up_score, left_score)

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
    The traceback stage.
    print a string of a path with arrows (i1j1 -> i2j2 -> i3j3). From the biggest number to 0, which is the starting point.
    If the three number around the biggest number are same, you can choose a random one to go.
    '''
    x, y = start_pos[0], start_pos[1]
    curr_point = matrix[x][y]
    loc = []
    loc.append((x,y))

    while curr_point != 0:
        if x-1 >= 0 and y-1 >= 0:
            potential = []
            x_up, y_up = x, y-1
            x_left, y_left = x-1, y
            x_dia, y_dia = x-1, y-1
            if curr_point - match == matrix[x_dia][y_dia]: # first check the diagnal number to see if they are either match or mismatch
                loc.append((x_dia, y_dia))
                x, y = x_dia, y_dia
            else: # after checking the diagnal number, take the max number from three directions (diagnal, up, left)
                if x_dia == 0 and y_dia == 0:
                    potential.append((matrix[x_up][y_up], x_up, y_up))
                    potential.append((matrix[x_dia][y_dia], x_dia, y_dia))
                    potential.append((matrix[x_left][y_left], x_left, y_left))
                else:
                    potential.append((matrix[x_dia][y_dia], x_dia, y_dia))
                    potential.append((matrix[x_up][y_up], x_up, y_up))
                    potential.append((matrix[x_left][y_left], x_left, y_left))
                target = max(potential, key = lambda x:x[0])
                curr_point = target[0]
                x, y = target[1], target[2]
                loc.append((x,y))
        else:
            break
    result = 'The path is: '
    for coord in loc:
        if coord == loc[-1]:
            result += str(coord)
        else:
            result += str(coord) + ' -> '
    print(result)
    return loc

def final_alignment(pos, matrix):
    '''
    Based on the traceback results, this func print out the final alignment results
    '''

    text1 = ''
    text2 = ''
    match_sign = []
    gap_count = 0

    # align the bases based on matrix score.
    curr_score_list = []
    for i in range(len(pos)):
        curr_score_list.append((matrix[pos[i][0]][pos[i][1]], pos[i][0], pos[i][1]))
    curr_score_list = curr_score_list + [(-200, -200, -200)]
    for i in range(len(curr_score_list)):
        if i < len(curr_score_list) - 1:
            if curr_score_list[i][0] == curr_score_list[i+1][0] + match: # match
                text1 += seq1[curr_score_list[i][1]-1]
                text2 += seq2[curr_score_list[i][2]-1]
                match_sign.append('*')
            elif curr_score_list[i][0] == curr_score_list[i+1][0] - match: # mismatch
                text1 += seq1[curr_score_list[i][1]-1]
                text2 += seq2[curr_score_list[i][2]-1]
                match_sign.append('.')
            else:   # gap
                if curr_score_list[i][2]-1 >= 0:
                    text1 += '-'
                    text2 += seq2[curr_score_list[i][2]-1]
                    match_sign.append(' ')
                    gap_count += 1



    text1 = text1[::-1]
    text2 = text2[::-1]
    match_sign = ''.join(match_sign)[::-1]

    # Adding not aligned bases to the string.
    seq1_len, seq2_len = len(seq1), len(seq2)
    beginning = list(pos[-1])
    ending = list(pos[0])
    temp_seq1_len, temp_seq2_len = seq1_len, seq2_len
    temp_seq1_len, temp_seq2_len = temp_seq1_len - len(pos) + gap_count, temp_seq2_len - len(pos)

    # If there is left over bases at the beginning of sequences.
    if beginning[0] != 0: # when there are bases left at the beginning of sequence 1
        b_diff = max(beginning[0], beginning[1])
        match_sign = ' '*b_diff + match_sign

        # deal with the bases at the beginning if sequence 1 has more leftovers.
        if beginning[0] == b_diff:
            temp1 = ''
            for i in range(0, beginning[0]):
                temp1 += seq1[i]
                temp_seq1_len -= 1
            text1 = temp1 + text1
        else:
            text1 = '-'*b_diff + text1

        # deal with the bases at the beginning if sequence 2 has more leftovers
        if beginning[1] == b_diff:
            temp2 = ''
            for i in range(0, beginning[1]):
                temp2 += seq2[i]
                temp_seq2_len -= 1
            text2 = temp2 + text2
        else:
            text2 = '-'*b_diff + text2

    # when there are no bases left at the beginning of sequence 1, but leftover at the beginning of sequence 2
    if beginning[0] == 0 and text2[0] != seq2[0] and text2[0] != '-':
        temp2 = ''
        for i in range(0, beginning[1]):
            temp2 += seq2[i]
            temp_seq2_len -= 1
        text2 = temp2 + text2
        text1 = '-'*beginning[1] + text1
        match_sign = ' '*beginning[1] + match_sign


    seq1_ending_count = 0 # count how many bases are added to the end of sequence 1
    seq2_ending_count = 0 # count how many bases are added to the end of sequence 2

    # adding left over bases to the end of sequence 1 and 2, if there is any.
    e_diff = max(temp_seq1_len, temp_seq2_len) # calculate which sequence has longer leftover bases in the end
    if temp_seq1_len == e_diff: # if sequence 1 has more, add more to the end of sequence 1
        for i in range(ending[0], seq1_len):
            text1 = text1 + seq1[i]
            seq1_ending_count += 1
    else: # if sequence 1 do not have more bases left than sequence 2, add whatever is left over for sequence 1, if any, and then add '-' to indicate gap in the end.
        for i in range(ending[0] + temp_seq1_len, seq1_len):
            text1 = text1 + seq1[i]
            seq1_ending_count += 1
        text1 = text1 + '-'*abs(e_diff - temp_seq1_len)
        seq1_ending_count += abs(e_diff - temp_seq1_len)

    if temp_seq2_len == e_diff: # if sequence 2 has more, add more to the end of sequence 2
        for i in range(ending[1], seq2_len):
            text2 = text2 + seq2[i]
            seq2_ending_count += 1
    else:# if sequence 2 do not have more bases left than sequence 1, add whatever is left over for sequence 2, if any, and then add '-' to indicate gap in the end.
        for i in range(ending[1] + temp_seq2_len, seq2_len):
            text2 = text2 + seq2[i]
            seq2_ending_count += 1
        text2 = text2 + '-'*abs(e_diff - temp_seq2_len)
        seq2_ending_count += abs(e_diff - temp_seq2_len)

    ending_adding = max(seq1_ending_count, seq2_ending_count) # calculate which seqence has more added bases in the end.

    # adding mismatch sign to the end of sequences.
    for i in range(len(text1)-ending_adding, len(text1)):
        if text1[i] != text2[i] and text1[i] != '-' and text2[i] != '-':
            match_sign += '.'
        else:
            match_sign += ' '

    # display the results
    print('text1:   ',text1)
    print('text2:   ',text2)
    print('         ',match_sign)

#end of your function(s)

if __name__ == '__main__':
    #my main
    score_matrix, start_pos = create_score_matrix(rows, cols)
    print_matrix(score_matrix)
    # find_path(score_matrix, start_pos)
    final_alignment(find_path(score_matrix, start_pos), score_matrix)

