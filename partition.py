from typing import Tuple, List

import numpy as np


def get_query_set(attributes: np.ndarray, usage: np.ndarray, strict: bool = True) -> List[int]:
    """
    Find set of queries that access attributes
    In strict mode, a given query must access all the attributes.
    When strict=False, the query may access any attribute in attributes

    :param attributes: attributes accessed by a query
    :param usage: attribute usage matrix
    :param strict: flag for specifying the type of constraint used
    :return: list of queries accessing attributes
    """
    def constraint(index: int):
        """
        Defines the constraint for strict and non-strict mode

        :param index: query being considered
        :return: True if constraint met, False otherwise
        """
        if strict:
            # usage has value 1 if a query is accessing that attribute
            # we need the query to access only those attributes
            # non-zero count in usage must be equal to number of attributes being considered
            return np.count_nonzero(usage[index]) == attributes.shape[0]
        else:
            # the query may access any attribute in attributes
            return True

    num_query = usage.shape[0]
    queries = []
    for i in range(num_query):
        flag = True
        for attr in attributes:
            # if ith query does not access attr, then ignore that query
            if usage[i, attr] != 1:
                flag = False
                break
        if flag and constraint(i):
            queries.append(i)

    return queries


def z_cost(x: int, ordering: np.ndarray, usage: np.ndarray, freq: np.ndarray, cost: np.ndarray) -> int:
    """
    Calculates cost for partition point x

    :param x: position to partition at
    :param ordering: attribute ordering after clustering
    :param usage: attribute usage matrix
    :param freq: query frequency matrix
    :param cost: query cost matrix
    :return: cost of partitioning at x
    """
    num_query, num_sites = freq.shape
    TA = ordering[:x + 1]
    BA = ordering[x + 1:]

    TQ = get_query_set(TA, usage)
    BQ = get_query_set(BA, usage)
    OQ = get_query_set(np.append(TA, BA), usage, False)

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


def get_partition_point(ordering: np.ndarray, usage: np.ndarray, freq: np.ndarray, cost: np.ndarray) -> int:
    """
    Find the partition point in the clustered attribute matrix (2-way partition)

    :param ordering: attribute order after clustering
    :param usage: attribute usage matrix
    :param freq: query frequency matrix
    :param cost: query cost matrix
    :return: index position of partition
    """
    num_attr = ordering.shape[0]

    costs = np.asarray([z_cost(x, ordering, usage, freq, cost) for x in range(1, num_attr - 1)])
    return np.argmax(costs) + 1


def get_shifted_point(ordering: np.ndarray, usage: np.ndarray, freq: np.ndarray, cost: np.ndarray) -> Tuple[int, int]:
    """
    Calculate the optimal partition point (2-way partition) after applying SHIFT operations

    :param ordering: attribute ordering after clustering
    :param usage: attribute usage matrix
    :param freq: query frequency matrix
    :param cost: query cost matrix
    :return: tuple with shift count and partition point
    """
    num_attr = ordering.shape[0]
    shift = 0
    point = 0
    best = -1_000_000

    # num_attr shifts can occur
    for i in range(num_attr):
        costs = np.asarray([z_cost(x, ordering, usage, freq, cost) for x in range(1, num_attr - 1)])
        if (m := np.amax(costs)) > best:
            best = m
            point = np.argmax(costs) + 1
            shift = i
        # the following operation "rolls" (shifts the elements) by -1 (to the left)
        # since ordering has one axis, it results in a simple left shift
        ordering = np.roll(ordering, -1)

    return shift, point
