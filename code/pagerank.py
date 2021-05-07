import time
import math
import os
import argparse
import numpy as np
from numpy.random import default_rng
from utils import *
from hypergraph import *

def read_data(dataset):
    print('dataset: {}'.format(dataset))
    H = HyperGraph()
    if dataset == 'dbpedia':
        H.read_file_pr('../data/dbpedia/dbpedia-writer_LCC.txt')
    elif dataset == 'kdd':
        H.read_file_pr('../data/kdd/dblp_kdd_LCC.txt')
    elif dataset == 'archive':
        H.read_file_pr('../data/archive/opsahl-collaboration_LCC.txt')
    elif dataset == 'highschool':
        H.read_cornell('highschool')
    elif dataset == 'ubuntu':
        H.read_cornell('ubuntu')
    H.dataset = dataset
    try:
        os.mkdir('../data/{}/tmp'.format(dataset))
    except:
        pass

    print('Input Hypergraph G:')
    H.origin_stat()
    return H

def sparsify_data(H, sparsify = 'site', c = 50, eps= 0.1, num_sites = 5):
    print('Parameter Settings:')
    print('c={}, eps={}, num_sites={}'.format(c, eps, num_sites))

    dataset = H.dataset
    if sparsify == 'site':
        print('\nSparsifying Local Hypergraphs:')
        H = H.sparsify_by_sites(num_sites, c, eps)
    else:
        print('Not sparsify')

    print('Union of Sparsifiers H in the Coordinator:')
    H.origin_stat()
    H.dataset = dataset
    return H


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', default='highschool', 
        choices = ['dbpedia', 'kdd', 'archive', 'highschool', 'ubuntu'])
    parser.add_argument('--num_sites', type=int, default=3)
    parser.add_argument('--c', type=int, default=2)
    args = parser.parse_args()

    c = args.c
    num_sites = args.num_sites
    dataset = args.dataset

    H = read_data(dataset)
    # print(H.ave_degree_weight())
    H1 = sparsify_data(H, sparsify='site', c = c, eps= 0.1, num_sites = num_sites)
    surffix = 'c{}_n{}'.format(c, num_sites)
    H.write(surffix)
    H1.write('sp' + surffix)

    rng = default_rng()
    numbers = rng.choice(H.n(), size=30, replace=False)

    print('\nRandomly sample 30 vertices to perform local clustering.')
    print('Vertex Indices: {}'.format(list(numbers)))
    print('Performing Clustering...')

    full_cond, full_t = [], []
    sp_cond, sp_t = [], []

    verbose = True if H.m() > 10000 else False
    for i,number in enumerate(numbers):
        if verbose:
            print('{}th V_init is {}'.format(i, number))
            print('Start work on Original Graph')
        cond, t = call_cs(dataset, '', number, c, num_sites, verbose)
        full_cond.append(cond)
        full_t.append(t)

        if verbose:
            print('Start work on Sparsified Graph')
        if len(H1.incident_edges[number]) == 0:
            print('Isolated point, conductance set as 0.')
            sp_cond.append(0)
            sp_t.append(0)
        else:
            cond, t = call_cs(dataset, 'sp', number, c, num_sites, verbose)
            sp_cond.append(cond)
            sp_t.append(t)

    print('Finish Clustering.')
    print('Clustering Quality:')
    print('Conductance Values in G: {}'.format(list(np.around(full_cond, decimals=4))))
    print('Conductance Values in H: {}'.format(list(np.around(sp_cond, decimals=4))))
    print('Clustering Time (second) in G: {}'.format(list(np.around(full_t, decimals=2))))
    print('Clustering Time (second) in H: {}'.format(list(np.around(sp_t, decimals=2))))
    print('---------')
    calculate_efficiency(H, H1, full_cond, sp_cond, full_t, sp_t)
