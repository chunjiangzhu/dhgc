import networkx as nx
import numpy as np
import subprocess
import seaborn as sns
import matplotlib.pyplot as plt
import os

PATH = os.getcwd()
# PageRank_exe = './Projects/HyperReplica/HyperReplica/bin/Release/HyperReplica.exe'
PageRank_exe = './HyperReplica.exe'

def check_connect(A):
    graph = nx.from_numpy_matrix(A)
    return nx.is_connected(graph)

def check_diag(A):
    diag_sum = np.diag(A).sum()
    return diag_sum == 0

def get_component_vertex(A):
    graph = nx.from_numpy_matrix(A)
    largest_cc = max(nx.connected_components(graph), key=len)
    return largest_cc

def calculate_efficiency(H, H1, full_cond, sp_cond, full_t, sp_t):
    full_cond, sp_cond = np.array(full_cond), np.array(sp_cond)
    full_t, sp_t = np.array(full_t), np.array(sp_t)
    ratio_comm2 =  ((H1.num_v_in_hyperedges())/H.num_v_in_hyperedges())
    ratio_abs_cond = (abs(full_cond - sp_cond) / full_cond).mean()
    print('Number of edges are {}, {}'.format(H.m(), H1.m()))
    print('Communication Costs are {}, {}'.format(H.num_v_in_hyperedges(), H1.num_v_in_hyperedges()))
    print('Ratio of Communication Cost is {}'.format( ratio_comm2))
    print('Conductance Quality abs: {}'.format(ratio_abs_cond))


def call_cs(dataset, sp, id, c, num_sites, verbose=True):
    surffix = 'c{}_n{}'.format(c, num_sites)
    path1 = '{}/{}_{}{}.txt'.format(PATH, dataset, sp, surffix)
    output1 = subprocess.check_output("mono {} {} {}".format(PageRank_exe, path1, id), shell=True)
    contents = (str(output1)).split('\\n')
    cond1 = float(contents[1].split(':')[1])
    time1 = float(contents[2].split(':')[1])
    if verbose:
        print('Conductance: {}, time {}'.format(cond1, time1))
    return cond1, time1

