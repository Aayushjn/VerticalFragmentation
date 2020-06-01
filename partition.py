from typing import Tuple

import numpy as np

from matrix import left_shift, up_shift


def z_cost(x: int, clustered: np.ndarray, usage: np.ndarray, freq: np.ndarray, cost: np.ndarray) -> int:
    """
    Calculates cost for partition point x

    :param x: position to partition at
    :param clustered: clustered attribute matrix
    :param usage: attribute usage matrix
    :param freq: query frequency matrix
    :param cost: query cost matrix
    :return: cost of partitioning at x
    """
    num_attr = clustered.shape[0]
    num_query, num_sites = freq.shape
    AQ, TQ, BQ, OQ = [], [], [], []

    for i in range(num_query):
        row = []
        for j in range(num_attr):
            if usage[i, j] == 1:
                row.append(j)
        AQ.append(row)

    start = num_attr - x
    for i in range(num_query):
        if AQ[i][1] <= start:
            TQ.append(i)
        elif AQ[i][0] > start:
            BQ.append(i)
        else:
            OQ.append(i)

    CTQ, CBQ, COQ = 0, 0, 0
    for i in range(len(TQ)):
        for j in range(num_sites):
            CTQ += (freq[TQ[i], j] * cost[TQ[i], j])
    for i in range(len(BQ)):
        for j in range(num_sites):
            CBQ += (freq[BQ[i], j] * cost[BQ[i], j])
    for i in range(len(OQ)):
        for j in range(num_sites):
            COQ += (freq[OQ[i], j] * cost[OQ[i], j])

    return (CTQ * CBQ) - (COQ ** 2)


def get_partition_point(clustered: np.ndarray, usage: np.ndarray, freq: np.ndarray, cost: np.ndarray) -> int:
    """
    Find the partition point in the clustered attribute matrix (2-way partition)

    :param clustered: clustered attribute matrix
    :param usage: attribute usage matrix
    :param freq: query frequency matrix
    :param cost: query cost matrix
    :return: index position of partition
    """
    num_attr = clustered.shape[0]

    costs = np.asarray([z_cost(x, clustered, usage, freq, cost) for x in range(1, num_attr - 1)])
    return np.argmax(costs) + 1


def get_shifted_point(clustered: np.ndarray, usage: np.ndarray, freq: np.ndarray, cost: np.ndarray) -> Tuple[int, int]:
    """
    Calculate the optimal partition point (2-way partition) after applying SHIFT operations

    :param clustered: clustered attribute matrix
    :param usage: attribute usage matrix
    :param freq: query frequency matrix
    :param cost: query cost matrix
    :return: tuple with shift count and partition point
    """
    num_attr = clustered.shape[0]
    shift = 0
    point = 0
    best = -1_000_000
    for i in range(num_attr):
        for j in range(num_attr - 1, 0, -1):
            z = z_cost(j, clustered, usage, freq, cost)
            if z > best:
                best = z
                point = j
                shift = i
        clustered = left_shift(up_shift(clustered))
        usage = left_shift(up_shift(usage))

    for i in range(shift):
        last = clustered[num_attr - 1]
        for j in range(num_attr - 1, 1, -1):
            clustered[j] = clustered[j - 1]
        clustered[0] = last

    return shift, point
