import h5py
import time
import pickle
import os
import subprocess
import numpy as np
import pandas as pd
from itertools import compress
import scipy.linalg
from scipy import sparse
from scipy.sparse import csr_matrix
from scipy.io import loadmat, savemat
from utils import *


class HyperGraph():
    def __init__(self):
        self.edges = None                           # List of list  (each edge include a list of vertex index)
        self.incident_edges = None                  # List of list  (each vertex include a list of edge index)
        self.weights = None                         # List
        self.ID = {}
        self.ID_rev = {}
        self.dataset = None

    def degree(self, v):
        # v is the index of vertex
        return len(self.incident_edges[v])

    def w_Degree(self, v):
        # v is the index of vertex
        sum = 0
        for e in self.incident_edges[v]:
            sum += self.weights[e]
        return  sum

    def n(self):
        return len(self.incident_edges)

    def m(self):
        return len(self.edges)

    def num_v_in_hyperedges(self):
        num = 0
        for edge in self.edges:
            num += len(edge)
        return num

    def ave_degree_weight(self):
        weight = 0
        for edge_list in self.incident_edges:
            for idx in edge_list:
                weight += self.weights[idx]
        return weight/self.n()

    def add_edge(self, edge, w=1):
        # here edge is a list of vertex index
        eid = self.m()
        self.edges.append(edge)
        for v in edge:
            while (v >= self.n()):
                self.incident_edges.append([])
            self.incident_edges[v].append(eid)
        self.weights.append(w)
        self.ID[edge] = eid
        # self.ID_rev[eid] = edge

    def read_file_pr(self, filename):
        # code from pagerank
        # each line include one edge (all the vextices indices belong to the edge) and associated weight.
        vertex_num = 0
        edges = []
        weights = []
        with open(filename, 'r') as fid:
            lines = fid.readlines()
            for line in lines:
                words = line.split()
                # i = 0
                edge = []
                for word in words[:-1]:
                    v = int(word)
                    edge.append(v)
                    if v >= vertex_num:
                        vertex_num = v + 1
                weights.append(float(words[-1]))
                # edge.sort()
                edges.append(edge)

        self.edges = edges
        self.weights = weights
        self.incident_edges = []
        for v in range(vertex_num):
            self.incident_edges.append([])
        for i in range(len(edges)):
            edge = edges[i]
            for v in edge:
                self.incident_edges[v].append(i)
            self.ID[tuple(edge)] = i
            self.ID_rev[i] = edge

    def write(self, surffix = ''):
        with open('../data/{}/tmp/{}_{}.txt'.format(self.dataset, self.dataset, surffix), 'w') as fid:
            for i in range(self.m()):
                line = ''
                edge = self.edges[i]
                for v in edge:
                    line += str(v)
                    line += ' '
                line += str(self.weights[i])
                line += '\n'
                fid.write(line)

    def read_cornell(self, dataset):
        if dataset == 'highschool':
            nverts_file = '../data/highschool/contact-high-school-nverts.txt'
            simplex_file = '../data/highschool/contact-high-school-simplices.txt'
        elif dataset == 'ubuntu':
            nverts_file = '../data/ubuntu/tags-ask-ubuntu-nverts.txt'
            simplex_file = '../data/ubuntu/tags-ask-ubuntu-simplices.txt'
        edge_num = np.array(pd.read_csv(nverts_file, header = None)).reshape(-1)
        edge_num_accumulate = np.cumsum(edge_num)

        simplex = np.array(pd.read_csv(simplex_file, header = None)).reshape(-1)
        edge_dic = {}
        edges = []
        weights = []
        for i in range(edge_num.shape[0]):
            if i == 0:
                start = 0
                end = edge_num_accumulate[i]
            else:
                start = end
                end = edge_num_accumulate[i]
            edge = list(simplex[start:end] - 1)
            edge.sort()
            edge = tuple(edge)
            if edge not in edge_dic:
                edge_dic[edge] = 1
            else:
                edge_dic[edge] += 1

        for key, value in edge_dic.items():
            edges.append(list(key))
            weights.append(value)

        vertex_num = simplex.max()
        self.edges = edges
        self.weights = weights
        self.incident_edges = []
        for v in range(vertex_num):
            self.incident_edges.append([])
        for i in range(len(edges)):
            edge = edges[i]
            for v in edge:
                self.incident_edges[v].append(i)
            self.ID[tuple(edge)] = i
            self.ID_rev[i] = edge

        self.to_component()
        adj = self.associated_graph()
        if not check_connect(adj):
            print('CHECK!! NOT CONNECTED.')

    def to_component(self):
        A = self.associated_graph()
        vertex = get_component_vertex(A)
        translate_dic = {}
        for i,v in enumerate(vertex):
            translate_dic[v] = i

        new_edges = []
        new_weights = []
        new_vertex_num = len(vertex)
        for i,edge in enumerate(self.edges):
            keep_edge = True
            new_edge = []
            for v in edge:
                if v not in vertex:
                    keep_edge = False
                    break
                else:
                    new_edge.append(translate_dic[v])
            if keep_edge:
                new_edges.append(new_edge)
                new_weights.append(self.weights[i])

        self.incident_edges = []
        self.ID = {}
        self.ID_rev = {}
        self.weights = new_weights
        self.edges = new_edges
        for v in range(new_vertex_num):
            self.incident_edges.append([])
        for i in range(len(new_edges)):
            edge = new_edges[i]
            for v in edge:
                self.incident_edges[v].append(i)
            self.ID[tuple(edge)] = i
            self.ID_rev[i] = edge

    def split_graph_sites(self, num_sites):
        file = '{}_{}.obj'.format(self.dataset, num_sites)
        path = '../data/{}/tmp'.format(self.dataset)
        files = os.listdir(path)
        if file in files:
            filehandler = open(os.path.join(path, file), 'rb')
            results = pickle.load(filehandler)
            return results

        sites = np.random.randint(num_sites, size=self.m())
        results = []
        for i in range(num_sites):
            flag = sites == i
            site_hypergraph = HyperGraph()
            site_hypergraph.edges = list(compress(self.edges, flag))
            site_hypergraph.weights = list(compress(self.weights, flag))
            site_hypergraph.incident_edges = []
            site_hypergraph.dataset = self.dataset
            for v in range(self.n()):
                site_hypergraph.incident_edges.append([])
            for j in range(len(site_hypergraph.edges)):
                edge = site_hypergraph.edges[j]
                for v in edge:
                    site_hypergraph.incident_edges[v].append(j)
                site_hypergraph.ID[tuple(edge)] = j
                site_hypergraph.ID_rev[j] = edge
            results.append(site_hypergraph)
        file_pi = open(os.path.join(path, file), 'wb')
        pickle.dump(results, file_pi)
        return results

    def sparsify_by_sites(self, num_sites, c, eps):
        graph_sites = self.split_graph_sites(num_sites)
        sparse_graph_sites = []
        total_time = 0
        total_num_edge = 0
        total_num_v = 0
        max_time = 0
        for i,h in enumerate(graph_sites):
            load = h.sparsify_stage1(num_sites, '{}'.format(i))
            if not load:
                output = subprocess.check_output("julia effective_resistance.jl {} {} {}".format(self.dataset, num_sites, i), shell=True)
                # print(output)
                time_re = float((str(output).split("\\n")[1]).split('=')[1])
            else:
                time_re = 0
            time_1 = time.time()
            h_sp = h.sparsify_stage2(c, eps, '{}'.format(i), num_sites)
            time_2 = time.time()
            time_sp = time_2 - time_1
            times = time_re + time_sp
            if times > max_time:
                max_time = times
            total_time += times
            print('Finished {}/{} site'.format(i+1, num_sites))
            total_num_edge += h_sp.m()
            total_num_v += h_sp.num_v_in_hyperedges()
            sparse_graph_sites.append(h_sp)
        H = self.merge_graphs(sparse_graph_sites)
        H.dataset = self.dataset
        print('Total time for Sparsifier is {}'.format(total_time))
        print('Max time for Sparsifier is {}'.format(max_time))
        return H

    def origin_stat(self):
        print('Total Number of Vertex is {}'.format(self.n()))
        print('Total Number of Hyperedge is {}'.format(self.m()))
        print('Total Number of vertex in Hyperedge is {}'.format(self.num_v_in_hyperedges()))

    def merge_graphs(self, graphs):
        H = HyperGraph()
        H.edges = [edge for g in graphs for edge in g.edges]
        H.weights = [weight for g in graphs for weight in g.weights]
        H.incident_edges = []
        H.ID = {}
        H.ID_rev = {}
        for v in range(self.n()):
            H.incident_edges.append([])
        for i in range(len(H.edges)):
            edge = H.edges[i]
            for v in edge:
                H.incident_edges[v].append(i)
            H.ID[tuple(edge)] = i
            H.ID_rev[i] = edge
        H.dataset = self.dataset
        return H

    def associated_graph(self):
        adj = np.zeros((self.n(), self.n()))
        for i in range(self.m()):
            edge = self.edges[i]
            weight = self.weights[i]
            for j in range(len(edge)):
                for k in range(j+1, len(edge)):
                    v1 = edge[j]
                    v2 = edge[k]
                    if v1 == v2:
                        print('One edge include same vertex: {}, {}'.format(v1, v2))
                    adj[v1, v2] += weight
                    adj[v2, v1] += weight
        return adj

    def indicent_matrix(self):
        I = np.zeros((self.m(), self.n()))
        for i, edge in enumerate(self.edges):
            for v in edge:
                I[i,v] = 1
        return I.transpose()

    def max_vertex_in_edges(self):
        r = 0
        for edge in self.edges:
            r = max(len(edge), r)
        return r

    def sparsify_stage1(self, num_sites, file_surfix):
        # In this stage, we save the associated graph and the pairs that need to calculoate the effective resistance
        # If the graph have been processed before, we directly return.
        files = os.listdir('../data/{}/tmp'.format(self.dataset))
        filename = 'pairs_n{}_{}.mat'.format(num_sites, file_surfix)
        self.pairs = self.pairs_for_er()
        if filename in files:
            return True

        A = self.associated_graph()
        A_sp = sparse.csr_matrix(A)
        savemat('../data/{}/tmp/sp_A_n{}_{}.mat'.format(self.dataset, num_sites, file_surfix), {'A_sp': A_sp})
        savemat('../data/{}/tmp/pairs_n{}_{}.mat'.format(self.dataset, num_sites, file_surfix), {'pairs': np.array(self.pairs)})
        return False

    def sparsify_stage2(self, c, eps, file_surfix, num_sites):
        # sparsify the graph according to the effective resistance calculated by julia
        # c is parameter need to be tuned
        # eps is the quality of the sparsifier
        content = h5py.File('../data/{}/tmp/rs_n{}_{}.mat'.format(self.dataset, num_sites, file_surfix))
        rs = content['rs']
        pairs = self.pairs

        rs_dic = {}
        for i in range(len(pairs)):
            key = tuple(pairs[i])
            try:
                rs_dic[key] = rs[i]
            except:
                pass
            else:
                pass

        # calculate r_e for hypergraph
        re_list = []
        for i in range(self.m()):
            max = -1e6
            edge = self.edges[i]
            for j in range(len(edge)):
                for k in range(j+1, len(edge)):
                    v1 = edge[j]
                    v2 = edge[k]
                    if v1 < v2:
                        key = tuple([v1, v2])
                    else:
                        key = tuple([v2, v1])
                    if rs_dic[key] > max:
                        max = rs_dic[key]
            re_list.append(max)

        # the formula includes r term, we count the r as a part of c to accelerate para tuning
        r = self.max_vertex_in_edges()
        L = eps **2 * c

        sampled_edges = []
        sampled_edges_weights = []
        for i in range(self.m()):
            edge = self.edges[i]
            re = re_list[i]
            pe = min(1, re/L)
            if np.random.uniform(0,1) < pe:
                sampled_edges.append(edge)
                sampled_edges_weights.append(self.weights[i] / pe)

        sparsified = HyperGraph()
        sparsified.edges = sampled_edges
        sparsified.weights = sampled_edges_weights
        sparsified.incident_edges = []
        for v in range(self.n()):
            sparsified.incident_edges.append([])
        for i in range(len(sparsified.edges)):
            edge = sparsified.edges[i]
            for v in edge:
                sparsified.incident_edges[v].append(i)
            sparsified.ID[tuple(edge)] = i
            sparsified.ID_rev[i] = edge

        return sparsified

    def pairs_for_er(self):
        # some of the pairs are replicated, we remove the duplicated pairs.
        pair_list = []
        for edge in self.edges:
            for i in range(len(edge)):
                for j in range(i+1, len(edge)):
                    if edge[i] < edge[j]:
                        pair = (edge[i], edge[j])
                    else:
                        pair = (edge[j], edge[i])
                    pair_list.append(pair)
        pair_list = list(set(pair_list))
        return pair_list










