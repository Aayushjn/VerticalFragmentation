from typing import List, Tuple


def get_attribute_affinity_matrix(usage: List[List[int]],
                                  freq: List[List[int]],
                                  cost: List[List[int]]) -> List[List[int]]:
    """
    Calculates attribute affinity matrix
    :param usage: attribute usage matrix
    :param freq: query frequency matrix
    :param cost: query cost matrix
    :return: attribute affinity matrix
    """
    num_attr = len(usage[0])
    num_sites = len(freq[0])
    num_query = len(freq)
    matrix = [[0 for _ in range(num_attr)] for _ in range(num_attr)]

    for i in range(num_attr):
        for j in range(num_attr):
            out_sum = 0
            for k in range(num_query):
                if usage[k][i] == 1 and usage[k][j] == 1:
                    in_sum = 0
                    for l in range(num_sites):
                        in_sum += (freq[k][l] * cost[k][l])
                    out_sum += in_sum
            matrix[i][j] = out_sum

    return matrix


def bond_energy_algorithm(affinity: List[List[int]]) -> Tuple[List[List[int]], List[int]]:
    """
    Apply bond energy algorithm and derive the clustered attribute matrix
    :param affinity: attribute affinity matrix
    :return: clustered attribute matrix and ordering of attributes
    """
    num_attr = len(affinity)
    matrix = [[0 for _ in range(num_attr)] for _ in range(num_attr)]

    # Fix 1st column in the clustered matrix and set current index
    order, index = [0], 1
    while index < num_attr:
        max_index = -1
        max_cont = -1_000_000
        for i in range(1, index):
            con = cont(order[i - 1], index, order[i], affinity)
            if con > max_cont:
                max_index = i
                max_cont = con

        # Set first attribute (-1 as ith element)
        con = cont(-1, index, order[0], affinity)
        if con > max_cont:
            max_index = 0
            max_cont = con

        # Set last attribute (-1 as jth element)
        con = cont(order[index - 1], index, -1, affinity)
        if con > max_cont:
            max_index = index

        # Set the final attribute
        if max_index == index:
            order.append(index)
        else:
            # Reorder attributes
            order.append(0)
            for j in range(index, max_index, -1):
                order[j] = order[j - 1]
            order[max_index] = index
        index += 1

    for i in range(num_attr):
        for j in range(num_attr):
            matrix[i][j] = affinity[order[i]][order[j]]

    # for i in range(1, num_attr):
    #     n_pos = len(order) + 1
    #     bond_strength = [0 for _ in range(num_attr)]
    #
    #     for p in range(n_pos):
    #         if p >= 1:
    #             bond_left = 2 * affinity[order[p - 1]][i]
    #         else:
    #             bond_left = 0
    #
    #         if p < (n_pos - 1):
    #             bond_right = 2 * affinity[order[p]][i]
    #         else:
    #             bond_right = 0
    #
    #         if 1 <= p < (n_pos - 1):
    #             bond_mid = 2 * affinity[order[p - 1]][order[p]]
    #         else:
    #             bond_mid = 0
    #
    #         bond_strength[p] = bond_left + bond_right - bond_mid
    #
    #     max_pos = -1
    #     max_val = -1000000
    #     for j in range(num_attr):
    #         if bond_strength[j] > max_val:
    #             max_val = bond_strength[j]
    #             max_pos = j
    #     P = [0 for _ in range(num_attr)]
    #     P[0:max_pos] = order[0:max_pos]
    #     P[max_pos] = i
    #     P[max_pos + 1:] = order[max_pos:]
    #
    #     order = P

    return matrix, order


def cont(i: int, k: int, j: int, affinity: List[List[int]]) -> int:
    """
    Calculate contribution of ordering affinity[i], affinity[k], affinity[j]
    :param i: ith column
    :param k: kth column
    :param j: jth column
    :param affinity: attribute affinity matrix
    :return: Contribution of order(i-k-j)
    """
    return (2 * bond(i, k, affinity)) + (2 * bond(k, i, affinity)) - (2 * bond(i, j, affinity))


def bond(x: int, y: int, affinity: List[List[int]]) -> int:
    """
    Calculate bond energy between affinity[x] and affinity[y]
    :param x: xth column
    :param y: yth column
    :param affinity: attribute affinity matrix
    :return: bond energy
    """
    # x or y may be -1 when determining contribution of ordering w.r.t index 0 or (num_attr - 1)
    if x == -1 or y == -1:
        return 0

    num_attr = len(affinity)
    bond_sum = 0
    for z in range(1, num_attr):
        bond_sum += (affinity[z][x] * affinity[z][y])
    return bond_sum
