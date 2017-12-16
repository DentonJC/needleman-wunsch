"""
https://en.wikipedia.org/wiki/Needleman–Wunsch_algorithm
"""

import numpy as np
from Bio.SubsMat.MatrixInfo import blosum62
from Bio import Entrez
Entrez.email = "your@mail.com"


def matrix(seq1, seq2, gap_penalty):
    """
    Create and fill score matrix
    """
    F = [[None] * (len(seq2)+1) for i in range(len(seq1)+1)]
    F[0][0] = 0
    for i in range(len(seq1)+1):
        F[i][0] = i * gap_penalty
    for i in range(len(seq2)+1):
        F[0][i] = i * gap_penalty

    for x in range(len(seq1)):
        for y in range(len(seq2)):
            i, j = x+1, y+1
            try:
                Match = F[i-1][j-1] + blosum62[(seq1[i-1], seq2[j-1])]
            except KeyError:
                Match = F[i-1][j-1] + blosum62[(seq2[j-1], seq1[i-1])]
            Delete = F[i-1][j] + gap_penalty
            Insert = F[i][j-1] + gap_penalty
            F[i][j] = max(Match, Delete, Insert)
    return np.array(F)


def NW(seq1, seq2, F):
    """
    Needleman–Wunsch algorithm
    """
    AlignmentA = ""
    AlignmentB = ""
    Score = 0
    i, j = len(seq1), len(seq2)
    while (i > 0 or j > 0):
        try:
            Match = F[i-1][j-1] + blosum62[(seq1[i-1], seq2[j-1])]
        except KeyError:
            Match = F[i-1][j-1] + blosum62[(seq2[j-1], seq1[i-1])]
        if i > 0 and j > 0 and F[i][j] == Match:
            Score += Match - F[i-1][j-1]
            AlignmentA = seq1[i-1] + AlignmentA
            AlignmentB = seq2[j-1] + AlignmentB
            i = i - 1
            j = j - 1
        elif(i > 0 and F[i][j] == F[i-1][j] + gap_penalty):
            Score += gap_penalty
            AlignmentA = seq1[i-1] + AlignmentA
            AlignmentB = "-" + AlignmentB
            i = i - 1
        else:
            Score += gap_penalty
            AlignmentA = "-" + AlignmentA
            AlignmentB = seq2[j-1] + AlignmentB
            j = j - 1
    return AlignmentA, AlignmentB, Score


if __name__ == "__main__":
    gap_penalty = -1
    id1 = '40886941'
    id2 = '34849618'

    handle = Entrez.efetch(db="protein", id=id1, rettype="gp", retmode="xml")
    seq1 = Entrez.read(handle)[0]["GBSeq_sequence"].upper()
    handle = Entrez.efetch(db="protein", id=id2, rettype="gp", retmode="xml")
    seq2 = Entrez.read(handle)[0]["GBSeq_sequence"].upper()

    F = matrix(seq1, seq2, gap_penalty)
    AlignmentA, AlignmentB, Score = NW(seq1, seq2, F)
    print(AlignmentA)
    print(AlignmentB)
    print("Score: ", Score)
