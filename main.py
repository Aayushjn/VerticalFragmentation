import sys

import numpy as np

from bea import get_attribute_affinity_matrix, bond_energy_algorithm


def main(file: str):
    try:
        with(open(file, 'r')) as f:
            num_attr = int(f.readline())
            num_query = int(f.readline())
            num_sites = int(f.readline())

            attr_usage_matrix = np.loadtxt(f, np.int, delimiter=' ', max_rows=num_query)
            query_freq_matrix = np.loadtxt(f, np.int, delimiter=' ', max_rows=num_query)
            query_cost_matrix = np.loadtxt(f, np.int, delimiter=' ', max_rows=num_query)
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

    # AQ, TQ, BQ, OQ = [], [], [], []
    # for i in range(num_query):
    #     row = []
    #     for j in range(num_attr):
    #         if attr_usage_matrix[i][j] == 1:
    #             row.append(j)
    #     AQ.append(row)
    #
    # for i in range(num_query):
    #     for j in range(len(AQ[i])):


if __name__ == '__main__':
    if len(sys.argv) == 1:
        print('Provide a .txt file as argument!')
        sys.exit(1)
    elif sys.argv[1].endswith('.txt'):
        main(sys.argv[1])
    else:
        print('Must be a .txt file!')
        sys.exit(2)
