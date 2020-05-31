import numpy as np


def up_shift(matrix: np.ndarray) -> np.ndarray:
    """
    Up shifts a square matrix by 1

    :param matrix: input matrix
    :return: up-shifted matrix
    """
    assert matrix.shape[0] == matrix.shape[1]

    return np.roll(matrix, -1, 0)


def left_shift(matrix: np.ndarray) -> np.ndarray:
    """
    Left shifts a square matrix by 1

    :param matrix: input matrix
    :return: left-shifted matrix
    """
    assert matrix.shape[0] == matrix.shape[1]

    return np.roll(matrix, -1, 1)
