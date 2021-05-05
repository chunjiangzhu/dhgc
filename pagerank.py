from hypergraph import *
import time
import math
import numpy as np
import os
from numpy.random import default_rng
from utils import *

def read_data(dataset):
    print('dataset: {}'.format(dataset))
    H = HyperGraph()
    if dataset == 'pr':
        H.read_file_pr('dbpedia-writer_LCC.txt')
    elif dataset == 'kdd':
        H.read_file_pr('dblp_kdd_LCC.txt')
    elif dataset == 'archive':
        H.read_file_pr('opsahl-collaboration_LCC.txt')
    elif dataset == 'highschool':
        H.read_cornell('highschool')
    elif dataset == 'ubuntu':
        H.read_cornell('ubuntu')
    H.dataset = dataset

    H.origin_stat()
    try:
        os.mkdir('./{}'.format(dataset))
    except:
        pass
    else:
        pass
    return H

def sparsify_data(H, sparsify = 'site', c = 50, eps= 0.1, num_sites = 5):
    print('c={}, eps={}, num_sites={}'.format(c, eps, num_sites))
    dataset = H.dataset
    if sparsify == 'site':
        H = H.sparsify_by_sites(num_sites, c, eps)
    else:
        print('Not sparsify')
    H.origin_stat()
    H.dataset = dataset
    return H

c = 2
num_sites = 3
dataset = 'highschool'

H = read_data(dataset)
print(H.ave_degree_weight())
H1 = sparsify_data(H, sparsify='site', c = c, eps= 0.1, num_sites = num_sites)
surffix = 'c{}_n{}'.format(c, num_sites)
H.write(surffix)
H1.write('sp' + surffix)

rng = default_rng()
numbers = rng.choice(H.n(), size=30, replace=False)
print('Vertex index are {}'.format(numbers))

full_cond, full_t = [], []
sp_cond, sp_t = [], []

verbose = True if H.m() > 10000 else False
for i,number in enumerate(numbers):
    if verbose:
        print('{}th V_init is {}'.format(i, number))

    print('Original Graph')
    cond, t = call_cs(dataset, '', number, c, num_sites, verbose)
    full_cond.append(cond)
    full_t.append(t)

    print('Sparsified Graph')
    if len(H1.incident_edges[number]) == 0:
        print('Isolated point, conductance set as 0.')
        sp_cond.append(0)
        sp_t.append(0)
    else:
        cond, t = call_cs(dataset, 'sp', number, c, num_sites, verbose)
        sp_cond.append(cond)
        sp_t.append(t)

print('Full Graph conductance: {}'.format(full_cond))
print('Sparse Graph conductance: {}'.format(sp_cond))
print('Full time: {}'.format(full_t))
print('Sparse time: {}'.format(sp_t))

calculate_efficiency(H, H1, full_cond, sp_cond, full_t, sp_t)
