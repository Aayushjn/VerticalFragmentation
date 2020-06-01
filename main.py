import sys

import numpy as np

from bea import get_attribute_affinity_matrix, bond_energy_algorithm
from partition import get_partition_point, get_shifted_point


def main(file: str):
    try:
        with(open(file, 'r')) as f:
            num_attr = int(f.readline())
            num_query = int(f.readline())
            num_sites = int(f.readline())

            attr_usage_matrix = np.loadtxt(f, np.int, delimiter=' ', max_rows=num_query)
            assert attr_usage_matrix.shape == (num_query, num_attr)

            query_freq_matrix = np.loadtxt(f, np.int, delimiter=' ', max_rows=num_query)
            assert query_freq_matrix.shape == (num_query, num_sites)

            query_cost_matrix = np.loadtxt(f, np.int, delimiter=' ', max_rows=num_query)
            assert query_cost_matrix.shape == (num_query, num_sites)
    except AssertionError:
        print(f'Usage matrix must be a ({num_query}x{num_attr}) matrix')
        print(f'Frequency matrix must be a ({num_query}x{num_sites}) matrix')
        print(f'Cost matrix must be a ({num_query}x{num_sites}) matrix')
        sys.exit(3)
    except FileNotFoundError:
        print(f'{file} does not exist!')
        sys.exit(4)

    attr_aff_matrix = get_attribute_affinity_matrix(attr_usage_matrix, query_freq_matrix, query_cost_matrix)
    print('Attribute Affinity Matrix')
    print(attr_aff_matrix)

    clustered_attr_matrix, ordering = bond_energy_algorithm(attr_aff_matrix)
    print('Clustered Attribute Matrix')
    print(clustered_attr_matrix)
    print(f'Attribute order: {ordering}')

    point = get_partition_point(ordering, attr_usage_matrix, query_freq_matrix, query_cost_matrix)

    print(f'Split point without SHIFT operation is {point}')

    TA = ordering[:point + 1]
    BA = ordering[point + 1:]

    print(f'TA = {TA}')
    print(f'BA = {BA}')

    shift, point = get_shifted_point(ordering, attr_usage_matrix, query_freq_matrix, query_cost_matrix)

    print(f'Split point with SHIFT = {shift} operation is {point}')

    # shift left by `-shift`
    ordering = np.roll(ordering, -shift)
    TA = ordering[:point + 1]
    BA = ordering[point + 1:]

    print(f'TA = {TA}')
    print(f'BA = {BA}')


if __name__ == '__main__':
    if len(sys.argv) == 1:
        print('Provide a .txt file as argument!')
        sys.exit(1)
    elif sys.argv[1].endswith('.txt'):
        main(sys.argv[1])
    else:
        print('Must be a .txt file!')
        sys.exit(2)
