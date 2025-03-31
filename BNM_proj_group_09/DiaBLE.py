#! /usr/bin/env python

"""
# -----------------------------------------------------------------------
# encoding: utf-8
# DIAMOnD.py
# Joerg Menche, Susan D. Ghiassian
# Last Modified: 2020-22-09
# Modified for DiaBLE support (custom universe_size)
#
# -----------------------------------------------------------------------
"""

import time
import networkx as nx
import numpy as np
import copy
import scipy.stats
from collections import defaultdict
import csv
import sys

# =============================================================================
def print_usage():
    print(' ')
    print('        usage: python3 DIAMOnD.py network_file seed_file n alpha(optional) outfile_name (optional)')
    print('        -----------------------------------------------------------------')
    print('        network_file : The edgelist must be provided as any delimiter-separated')
    print('                       table. Make sure the delimiter does not exist in gene IDs')
    print('                       and is consistent across the file.')
    print('                       The first two columns of the table will be')
    print('                       interpreted as an interaction gene1 <==> gene2')
    print('        seed_file    : table containing the seed genes (if table contains')
    print('                       more than one column they must be tab-separated;')
    print('                       the first column will be used only)')
    print('        n            : desired number of DIAMOnD genes, 200 is a reasonable')
    print('                       starting point.')
    print('        alpha        : an integer representing weight of the seeds, default')
    print('                       value is set to 1')
    print('        outfile_name : results will be saved under this file name')
    print('                       by default the outfile_name is set to "first_n_added_nodes_weight_alpha.txt"')
    print(' ')

# =============================================================================
def check_input_style(input_list):
    try:
        network_edgelist_file = input_list[1]
        seeds_file = input_list[2]
        max_number_of_added_nodes = int(input_list[3])
    except:
        print_usage()
        sys.exit(0)
        return

    alpha = 1
    outfile_name = 'first_%d_added_nodes_weight_%d.txt' % (max_number_of_added_nodes, alpha)

    if len(input_list) == 5:
        try:
            alpha = int(input_list[4])
            outfile_name = 'first_%d_added_weight_%d.txt' % (max_number_of_added_nodes, alpha)
        except:
            outfile_name = input_list[4]

    if len(input_list) == 6:
        try:
            alpha = int(input_list[4])
            outfile_name = input_list[5]
        except:
            print_usage()
            sys.exit(0)
            return

    return network_edgelist_file, seeds_file, max_number_of_added_nodes, alpha, outfile_name

# =============================================================================
def read_input(network_file, seed_file):
    """
    Reads the network and the list of seed genes from external files.
    * The edgelist must be provided as a delimiter-separated table.
    * The seed genes must be provided as a table. The first column is used only.
    * Lines that start with '#' will be ignored in both files.
    """
    sniffer = csv.Sniffer()
    line_delimiter = None
    for line in open(network_file, 'r'):
        if line.startswith('#'):
            continue
        else:
            dialect = sniffer.sniff(line)
            line_delimiter = dialect.delimiter
            break

    if line_delimiter is None:
        print("network_file format not correct")
        sys.exit(0)

    # read the network:
    G = nx.Graph()
    for line in open(network_file, 'r'):
        if line.startswith('#'):
            continue
        line_data = line.strip().split(line_delimiter)
        node1, node2 = line_data[0], line_data[1]
        G.add_edge(node1, node2)

    # read the seed genes:
    seed_genes = set()
    for line in open(seed_file, 'r'):
        if line.startswith('#'):
            continue
        line_data = line.strip().split('\t')
        seed_gene = line_data[0]
        seed_genes.add(seed_gene)

    return G, seed_genes

# ================================================================================
def compute_all_gamma_ln(N):
    """
    precomputes log gamma(i) for i up to N
    """
    gamma_ln = {}
    for i in range(1, N + 1):
        gamma_ln[i] = scipy.special.gammaln(i)
    return gamma_ln

# =============================================================================
def logchoose(n, k, gamma_ln):
    if n - k + 1 <= 0:
        return scipy.infty
    lgn1 = gamma_ln[n + 1]
    lgk1 = gamma_ln[k + 1]
    lgnk1 = gamma_ln[n - k + 1]
    return lgn1 - (lgnk1 + lgk1)

# =============================================================================
def gauss_hypergeom(x, r, b, n, gamma_ln):
    """
    hypergeometric pmf = C(r, x)*C(b, n-x) / C(r+b, n)
    """
    return np.exp( logchoose(r, x, gamma_ln)
                 + logchoose(b, n - x, gamma_ln)
                 - logchoose(r + b, n, gamma_ln) )

# =============================================================================
### MODIFIED pvalue(...) to accept universe_size instead of N
def pvalue(kb, k, universe_size, s, gamma_ln):
    """
    Computes the p-value for a node that has kb out of k links to seeds,
    given s seeds in a 'universe' of size universe_size.

    pval = sum_{n=kb to k} [ C(s, n)*C(universe_size - s, k-n) / C(universe_size, k) ]
    """
    p = 0.0
    for n in range(kb, k + 1):
        if n > s:
            break
        prob = gauss_hypergeom(n, s, universe_size - s, k, gamma_ln)
        p += prob
    return min(1.0, p)

def get_neighbors_and_degrees(G):
    neighbors = {}
    all_degrees = {}
    for node in G.nodes():
        neighbors[node] = set(G.neighbors(node))
        all_degrees[node] = G.degree(node)
    return neighbors, all_degrees

def reduce_not_in_cluster_nodes(all_degrees, neighbors, G, not_in_cluster, cluster_nodes, alpha):
    reduced_not_in_cluster = {}
    kb2k = defaultdict(dict)
    for node in not_in_cluster:
        k = all_degrees[node]
        kb = 0
        for neigh in neighbors[node]:
            if neigh in cluster_nodes:
                kb += 1

        # Weighted edges if alpha > 1
        k  += (alpha - 1) * kb
        kb += (alpha - 1) * kb
        kb2k[kb][k] = node

    # We'll pick for each kb the node with the smallest k,
    # then pick for each k the node with the largest kb:
    k2kb = defaultdict(dict)
    for _kb, k2node in kb2k.items():
        min_k = min(k2node.keys())
        node = k2node[min_k]
        k2kb[min_k][_kb] = node

    for k_val, kb2node in k2kb.items():
        max_kb = max(kb2node.keys())
        node = kb2node[max_kb]
        reduced_not_in_cluster[node] = (max_kb, k_val)

    return reduced_not_in_cluster

# ======================================================================================
#   C O R E    A L G O R I T H M
# ======================================================================================
### MODIFIED to accept universe_size
def diamond_iteration_of_first_X_nodes(G, S, X, alpha, universe_size=None):
    """
    Parameters:
    ----------
    G : graph
    S : seeds
    X : # of iterations (i.e. how many nodes to add)
    alpha : seed weight
    universe_size : custom universe size for hypergeometric test
                    (if None, defaults to G.number_of_nodes())

    Returns:
    --------
    added_nodes : list of (node_name, k, kb, p_value) in the order they are added.
    """

    # If user didn't specify, revert to the old approach (N = graph size)
    if universe_size is None:  ### MODIFIED
        universe_size = G.number_of_nodes()

    # We'll still keep track of N_graph for alpha weighting logic
    N_graph = G.number_of_nodes()

    # Setting up neighbors and degrees
    neighbors, all_degrees = get_neighbors_and_degrees(G)

    # cluster_nodes = seeds
    cluster_nodes = set(S)
    not_in_cluster = set()
    s0 = len(cluster_nodes)

    # Weighted approach:
    s0       += (alpha - 1) * s0
    N_graph  += (alpha - 1) * s0

    # Precompute gamma ln for the 'universe_size'
    gamma_ln = compute_all_gamma_ln(int(universe_size) + 1)

    # Initialize not_in_cluster
    for node in cluster_nodes:
        not_in_cluster |= neighbors[node]
    not_in_cluster -= cluster_nodes

    added_nodes = []
    all_p = {}

    while len(added_nodes) < X:
        reduced_dict = reduce_not_in_cluster_nodes(
            all_degrees, neighbors, G,
            not_in_cluster,
            cluster_nodes,
            alpha
        )

        pmin = 10
        next_node = None
        info_dict = {}

        for node, (kb, k) in reduced_dict.items():
            # We'll look up or compute the p-value
            if (k, kb, s0) not in all_p:
                # old code was p = pvalue(kb, k, N_graph, s0, gamma_ln)
                p_temp = pvalue(kb, k, universe_size, s0, gamma_ln)  ### MODIFIED
                all_p[(k, kb, s0)] = p_temp

            p = all_p[(k, kb, s0)]
            info_dict[node] = (k, kb, p)

            if p < pmin:
                pmin = p
                next_node = node

        if next_node is None:
            # No more nodes can be added
            break

        best_k, best_kb, best_p = info_dict[next_node]
        added_nodes.append((next_node, best_k, best_kb, best_p))

        # Now we update cluster_nodes
        cluster_nodes.add(next_node)
        s0 = len(cluster_nodes)

        # Expand not_in_cluster
        not_in_cluster |= (neighbors[next_node] - cluster_nodes)
        if next_node in not_in_cluster:
            not_in_cluster.remove(next_node)

    return added_nodes

# ===========================================================================
#
#   M A I N    D I A M O n D    A L G O R I T H M
#
# ===========================================================================
### MODIFIED to accept universe_size as well
def DIAMOnD(G_original, seed_genes, max_number_of_added_nodes, alpha,
            outfile=None, universe_size=None):
    """
    Runs the DIAMOnD (or DiaBLE) algorithm with the possibility of
    specifying a custom universe_size.

    Input:
    ------
     - G_original : networkx Graph
     - seed_genes : set of seed nodes
     - max_number_of_added_nodes : integer
     - alpha : seeds weight
     - outfile : if provided, results are saved here
     - universe_size : if provided, uses this # in hypergeometric test
                       instead of G_original.number_of_nodes()

    Returns:
    --------
      - added_nodes : list of 4-tuples: (node, k, kb, p)
    """

    # 1. remove seed genes that are not in the network
    all_genes_in_network = set(G_original.nodes())
    seed_genes = set(seed_genes)
    disease_genes = seed_genes & all_genes_in_network

    missing = len(seed_genes - all_genes_in_network)
    if missing > 0:
        print(f"DIAMOnD(): ignoring {missing} of {len(seed_genes)} seeds not in the network")

    # 2. run the iterative algorithm
    added_nodes = diamond_iteration_of_first_X_nodes(
        G_original,
        disease_genes,
        max_number_of_added_nodes,
        alpha,
        universe_size=universe_size
    )

    # 3. save results
    if outfile is not None:
        with open(outfile, 'w') as fout:
            fout.write('\t'.join(['#rank', 'DIAMOnD_node', 'p_hyper']) + '\n')
            for rank, (node, k, kb, p) in enumerate(added_nodes, start=1):
                fout.write('\t'.join([str(rank), str(node), str(p)]) + '\n')

    return added_nodes

# ===========================================================================
#
# "Hey Ho, Let's go!" -- The Ramones (1976)
#
# ===========================================================================

if __name__ == '__main__':
    input_list = sys.argv
    network_edgelist_file, seeds_file, max_number_of_added_nodes, alpha, outfile_name = check_input_style(input_list)
    G_original, seed_genes = read_input(network_edgelist_file, seeds_file)

    # run DIAMOnD with default 'universe_size=None'
    added_nodes = DIAMOnD(
        G_original,
        seed_genes,
        max_number_of_added_nodes,
        alpha,
        outfile=outfile_name,
        universe_size=None  # same as old DIAMOnD
    )

    print("\n results have been saved to '%s' \n" % outfile_name)
