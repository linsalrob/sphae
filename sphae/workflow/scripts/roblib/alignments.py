__author__ = 'Rob Edwards'

def ed_neighborsof(i,j, scorematrix):
    neighbors = []

    if i - 1 >= 0 and j - 1 >= 0:
        neighbors.append(scorematrix[i-1][j-1])
        neighbors.append(scorematrix[i-1][j])
        neighbors.append(scorematrix[i][j-1])
    elif i - 1 >= 0:
        neighbors.append(scorematrix[i-1][j])
        #neighbors.append(j+1)
    elif j - 1 >= 0:
        neighbors.append(scorematrix[i][j-1])
        #neighbors.append(i+1)
    else:
        neighbors.append(0)

    return neighbors

def ed_penalty(acid1, acid2):
    return 1


def edit_distance(seq1, seq2):

    scorematrix = [[0 for s in seq1] for t in seq2]

    for i in range(0, len(seq2)):
        for j in range(0, len(seq1)):
            if seq2[i] == seq1[j]:
                if i > 0 and j > 0:
                    scorematrix[i][j] = scorematrix[i-1][j-1]
                else:
                    scorematrix[i][j]= min(ed_neighborsof(i,j,scorematrix))
            else:
                scorematrix[i][j] = min(ed_neighborsof(i,j,scorematrix)) + ed_penalty(seq2[i], seq1[j])

    return scorematrix[len(seq2) - 1][len(seq1) - 1]

