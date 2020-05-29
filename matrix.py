import numpy as np


def up_shift(matrix: np.ndarray) -> np.ndarray:
    """
    Up shifts a square matrix by 1
    :param matrix: input matrix
    :return: up-shifted matrix
    """
    assert matrix.shape[0] == matrix.shape[1]

    rows = matrix.shape[0]
    first_row = [matrix[0, i] for i in range(rows)]

    for i in range(1, rows):
        for j in range(rows):
            matrix[i - 1, j] = matrix[i, j]
    for j in range(rows):
        matrix[rows - 1, j] = first_row[j]

    return matrix


def left_shift(matrix: np.ndarray) -> np.ndarray:
    """
    Left shifts a square matrix by 1
    :param matrix: input matrix
    :return: left-shifted matrix
    """
    assert matrix.shape[0] == matrix.shape[1]

    rows = matrix.shape[0]
    first_col = [matrix[i, 0] for i in range(rows)]

    for i in range(rows):
        for j in range(1, rows):
            matrix[i, j - 1] = matrix[i, j]
    for i in range(rows):
        matrix[i, rows - 1] = first_col[i]

    return matrix
