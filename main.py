from sys import argv, exit

from bea import get_attribute_affinity_matrix, bond_energy_algorithm
from matrix import print_matrix


def main(file: str):
    try:
        with(open(file, 'r')) as f:
            num_attr = int(f.readline())
            num_query = int(f.readline())
            num_sites = int(f.readline())

            attr_usage_matrix = [list(map(int, f.readline().split(' '))) for _ in range(num_query)]
            query_freq_matrix = [list(map(int, f.readline().split(' '))) for _ in range(num_query)]
            query_cost_matrix = [list(map(int, f.readline().split(' '))) for _ in range(num_query)]
    except FileNotFoundError:
        print(f'{file} does not exist!')
        exit(4)

    attr_aff_matrix = get_attribute_affinity_matrix(attr_usage_matrix, query_freq_matrix, query_cost_matrix)
    print('Attribute Affinity Matrix')
    print_matrix(attr_aff_matrix)

    clustered_attr_matrix, ordering = bond_energy_algorithm(attr_aff_matrix)
    print('Clustered Attribute Matrix')
    print_matrix(clustered_attr_matrix)
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
    #     if AQ[i][1] <= num_attr - 2:
    #         TQ.append(i)
    #     elif AQ[i][0] > num_attr - 2:
    #         BQ.append(i)
    #     else:
    #         OQ.append(i)
    #
    # print(AQ)
    # print(TQ)
    # print(BQ)
    # print(OQ)


if __name__ == '__main__':
    if len(argv) == 1:
        print('Provide a .txt file as argument!')
        exit(1)
    elif argv[1].endswith('.txt'):
        main(argv[1])
    else:
        print('Must be a .txt file!')
        exit(2)
