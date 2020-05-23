from typing import List


def print_matrix(matrix: List[List[int]]):
    """
    Prints a matrix row by row
    :param matrix: matrix to be printed
    """
    rows = len(matrix)
    cols = len(matrix[0])
    for i in range(rows):
        for j in range(cols):
            print(matrix[i][j], end=' ')
        print()
    print()


def up_shift(matrix: List[List[int]]) -> List[List[int]]:
    """
    Up shifts a square matrix by 1
    :param matrix: input matrix
    :return: up-shifted matrix
    """
    assert len(matrix) == len(matrix[0])

    rows = len(matrix)
    first_row = [matrix[0][i] for i in range(rows)]

    for i in range(1, rows):
        for j in range(rows):
            matrix[i - 1][j] = matrix[i][j]
    for j in range(rows):
        matrix[rows - 1][j] = first_row[j]

    return matrix


def left_shift(matrix: List[List[int]]) -> List[List[int]]:
    """
    Left shifts a square matrix by 1
    :param matrix: input matrix
    :return: left-shifted matrix
    """
    assert len(matrix) == len(matrix[0])

    rows = len(matrix)
    first_col = [matrix[i][0] for i in range(rows)]

    for i in range(rows):
        for j in range(1, rows):
            matrix[i][j - 1] = matrix[i][j]
    for i in range(rows):
        matrix[i][rows - 1] = first_col[i]

    return matrix
