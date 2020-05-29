from typing import Tuple

import numpy as np


def get_attribute_affinity_matrix(usage: np.ndarray,
                                  freq: np.ndarray,
                                  cost: np.ndarray) -> np.ndarray:
    """
    Calculates attribute affinity matrix
    :param usage: attribute usage matrix
    :param freq: query frequency matrix
    :param cost: query cost matrix
    :return: attribute affinity matrix
    """
    num_attr = usage.shape[1]
    matrix = np.zeros((num_attr, num_attr), np.int)

    for i in range(num_attr):
        for j in range(num_attr):
            matrix[i, j] = aff(i, j, usage, freq, cost)

    return matrix


def aff(i: int, j: int, usage: np.ndarray, freq: np.ndarray, cost: np.ndarray) -> int:
    """
    Calculates the number of times attributes i and j are accessed together, taking all sites into consideration
    :param i: ith attribute
    :param j: jth attribute
    :param usage: attribute usage matrix
    :param freq: query frequency matrix
    :param cost: query cost matrix
    :return: number of times attributes i and j were accessed together
    """
    num_query, num_sites = freq.shape
    out_sum = 0
    for k in range(num_query):
        if usage[k, i] == 1 and usage[k, j] == 1:
            in_sum = 0
            for l in range(num_sites):
                in_sum += (freq[k, l] * cost[k, l])
            out_sum += in_sum
    return out_sum


def bond_energy_algorithm(affinity: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Apply bond energy algorithm and derive the clustered attribute matrix
    :param affinity: attribute affinity matrix
    :return: clustered attribute matrix and ordering of attributes
    """
    num_attr = affinity.shape[0]
    clustered_matrix = np.zeros((num_attr, num_attr), np.int)

    order = np.zeros(num_attr, np.int)
    order[1] = 1
    index = 2
    while index < num_attr:
        max_index = -1
        max_cont = -1_000_000
        for i in range(index):
            con = cont(order[i - 1], index, order[i], affinity)
            if con > max_cont:
                max_index = i
                max_cont = con

        con = cont(order[index - 1], index, index + 1, affinity)
        if con > max_cont:
            max_index = index

        for j in range(index, max_index, -1):
            order[j] = order[j - 1]
        order[max_index] = index
        index += 1

    # Order columns
    for i in range(num_attr):
        clustered_matrix[:, i] = affinity[:, order[i]]

    # Order rows
    temp = clustered_matrix.copy()
    for i in range(num_attr):
        clustered_matrix[i, :] = temp[order[i], :]

    return clustered_matrix, order


def cont(i: int, k: int, j: int, affinity: np.ndarray) -> int:
    """
    Calculate contribution of ordering affinity[i], affinity[k], affinity[j]
    :param i: ith column
    :param k: kth column
    :param j: jth column
    :param affinity: attribute affinity matrix
    :return: Contribution of order(i-k-j)
    """
    return 2 * (bond(i, k, affinity) + bond(k, j, affinity) - bond(i, j, affinity))


def bond(x: int, y: int, affinity: np.ndarray) -> int:
    """
    Calculate bond energy between affinity[x] and affinity[y]
    :param x: xth column
    :param y: yth column
    :param affinity: attribute affinity matrix
    :return: bond energy
    """
    # x or y may be -1 when determining contribution of ordering w.r.t index 0 or (num_attr - 1)
    num_attr = affinity.shape[0]
    if x < 0 or y < 0 or x >= num_attr or y >= num_attr:
        return 0

    bond_sum = 0
    for z in range(num_attr):
        bond_sum += (affinity[z, x] * affinity[z, y])
    return bond_sum
