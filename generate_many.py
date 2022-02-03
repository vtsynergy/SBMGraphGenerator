import argparse
from collections import namedtuple
import getpass
import io
import os
import random

import graph_tool.all as gt
import numpy as np
import pandas as pd  # for writing output graph TSV files
import paramiko
import scipy.stats as stats


def analyze_generated_graph(graph: gt.Graph, block_membership: np.array, num_blocks: int):
    # compute and report basic statistics on the generated graph
    bg, bb, vcount, ecount, avp, aep = gt.condensation_graph(graph, block_membership, self_loops=True)
    edge_count_between_blocks = np.zeros((num_blocks, num_blocks))
    for e in bg.edges():
        edge_count_between_blocks[bg.vertex_index[e.source()], bg.vertex_index[e.target()]] = ecount[e]
    num_within_block_edges = sum(edge_count_between_blocks.diagonal())
    num_between_block_edges = graph.num_edges() - num_within_block_edges
    # print count statistics
    print('Number of nodes: {} Number of edges: {}'.format(graph.num_vertices(), graph.num_edges()))
    degrees = graph.get_total_degrees(np.arange(graph.num_vertices()))
    print('Vertex degrees: [{},{},{}]'.format(np.min(degrees), np.mean(degrees), np.max(degrees)))
    unique_degrees, counts = np.unique(degrees, return_counts=True)
    print("degrees: {}\ncounts: {}".format(unique_degrees[:20], counts[:20]))
    print('Avg. Number of nodes per block: {}'.format(graph.num_vertices() / num_blocks))
    print('# Within-block edges / # Between-blocks edges: {}'.format(num_within_block_edges / num_between_block_edges))
    # save_graph(graph, true_partition, utility, file_name)
# End of analyze_generated_graph()


def generate_block_membership(_num_blocks: int, num_vertices: int, block_size_heterogeneity: float) -> [int, np.array]:
    block_distribution = np.random.dirichlet(np.ones(_num_blocks) * 10 / block_size_heterogeneity, 1)[0]
    # draw block membership for each node
    block_membership_vector = np.where(np.random.multinomial(n=1, size=num_vertices, pvals=block_distribution))[1]
    # renumber this in case some blocks don't have any elements
    blocks, counts = np.unique(block_membership_vector, return_counts=True)
    block_mapping = {value: index for index, value in enumerate(blocks)}
    block_membership_vector = np.asarray([
        block_mapping[block_membership_vector[i]] for i in
        range(block_membership_vector.size)
    ])
    num_blocks = blocks.size
    return num_blocks, block_membership_vector
# End of generate_block_membership()


def generate_filename(num_vertices: int, num_blocks: int, max_degree: float, power_law_exponent: float,
                      block_overlap: float, block_size_variation: float, density: float,
                      args: argparse.Namespace) -> str:
    tag = "difficulty_{0}_{1}_{2}_{3}".format(num_blocks, max_degree, power_law_exponent, density)
    file_name = "{4}/{3}/{0}Overlap_{1}BlockSizeVar/{3}_{0}Overlap_{1}BlockSizeVar_{2}_nodes".format(
        block_overlap, block_size_variation, num_vertices, tag, args.directory)
    if args.remote:
        utility.exec_command("mkdir -p {}".format(os.path.dirname(file_name)))
    else:
        utility.makedirs(os.path.dirname(file_name), exist_ok=True)
    return file_name
# End of generate_filename()


def fill_blockmodel(num_vertices: int, bm_size: int, _overlap: float, block_membership_vector: np.array) -> np.array:
    # set the within-block and between-block edge strength accordingly
    def inter_block_strength(a, b):
        if a == b:  # within block interaction strength
            return 1
        else:  # between block interaction strength
            avg_within_block_nodes = float(num_vertices) / bm_size
            avg_between_block_nodes = num_vertices - avg_within_block_nodes
            return avg_within_block_nodes / avg_between_block_nodes / _overlap

    blocks, counts = np.unique(block_membership_vector, return_counts=True)
    blockmodel = np.zeros((bm_size, bm_size), dtype=np.float32)
    for row in range(bm_size):
        for col in range(bm_size):
            strength = inter_block_strength(row, col)
            value = strength * counts[row] * counts[col]
            blockmodel[row, col] = value
    return blockmodel
# End of fill_blockmodel()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--numvertices", type=int, default=200, help="Number of vertices in graph")
    # parser.add_argument("-c", "--communities", type=int, default=-1,
    #                     help="Number of communities in graph. If -1, uses communityexponent argument instead.")
    parser.add_argument("-m", "--communityexponent", type=float, default=0.35,
                        help="Number of communities = n^m. Only used if communities is -1.")
    # parser.add_argument("-i", "--mindegree", type=int, default=10, help="Min vertex degree = min(i, n / 4c)")
    parser.add_argument("-a", "--maxdegree", type=float, default=0.05, help="Max vertex degree = a * n")
    parser.add_argument("-o", "--overlap", type=float, default=5.0, help="5 = low, 2 (or 1.7) = high")
    parser.add_argument("-s", "--blocksizevariation", type=float, default=1.0, help="1 = low, 5 = high")
    parser.add_argument("-e", "--powerlawexponent", type=float, default=-2.1, help="The power law exponent")
    parser.add_argument("-d", "--density", type=float, default=1.0, help="1 - d fraction of edges will be removed")
    parser.add_argument("--directory", type=str, default="./", help="The directory in which to store the graphs")
    parser.add_argument("-r", "--remote", type=str, default="",
                        help="""If provided, will save files on remote system. Needs username and password to be
                             manually entered.""")
    return parser.parse_args()


def save_graph(graph, true_partition, utility, filename):
    """Saves graph using given utility is os, uses os.open() to write new file. Else, uses paramiko's
    sftp.open() to write new file. Vertices are indexed starting from 1.
    """
    # output the graph and truth partition to TSV files with standard format
    # store edge list
    edge_list = np.zeros((graph.num_edges(), 3), dtype=int)
    # populate the edge list.
    counter = 0
    # iterate through all edges (edge list access is not available in all versions of graph-tool)
    for e in graph.edges():
        edge_list[counter, 0] = int(e.source()) + 1  # nodes are indexed starting at 1 in the standard format
        edge_list[counter, 1] = int(e.target()) + 1  # nodes are indexed starting at 1 in the standard format
        edge_list[counter, 2] = 1  # all edges are weighted equally at 1 in this generator
        counter += 1
    # write graph TSV file using pandas DataFrame
    df_graph = pd.DataFrame(edge_list)
    graph_data = df_graph.to_csv(sep='\t', header=False, index=False)
    # write truth partition TSV file using pandas DataFrame; nodes and blocks are indexed starting at 1
    # in the standard format
    df_partition = pd.DataFrame(np.column_stack((np.arange(graph.num_vertices()) + 1, true_partition + 1)))
    partition_data = df_partition.to_csv(sep='\t', header=False, index=False)
    # print("utility type: ", type(utility))
    if isinstance(utility, paramiko.client.SSHClient): # str(type(utility)) == "paramiko.client.SSHClient":
        sftp = utility.open_sftp()
        with sftp.open("{}.tsv".format(filename), 'w') as file:
            file.write(graph_data)
        with sftp.open("{}_truePartition.tsv".format(filename), 'w') as file:
            file.write(partition_data)
    else:
        print(utility)
        print(filename)
        with io.open("{}.tsv".format(filename), 'w') as file:
            file.write(graph_data)
        with io.open("{}_truePartition.tsv".format(filename), 'w') as file:
            file.write(partition_data)
# End of save_graph()


# Generate the graph according to the blockmodel and parameters
def generate(args, utility):
    N_adjusted = int(args.numvertices * 1.13)
    M = args.communityexponent
    min_degree = 1  # args.mindegree
    max_degree = int(args.maxdegree * N_adjusted)  # A = args.maxdegree
    ratio_within_over_between = args.overlap
    block_size_heterogeneity = args.blocksizevariation
    powerlaw_exponent = args.powerlawexponent
    density = args.density
    num_blocks = int(N_adjusted ** M)  # number of blocks grows sub-linearly with number of nodes. Exponent (M) is a parameter.
    print('Number of blocks: {}'.format(num_blocks))
    # file_name = generate_filename(N, num_blocks, args.maxdegree, powerlaw_exponent, args.overlap,
    #                               args.blocksizevariation, density, args)

    # define discrete power law distribution
    def discrete_power_law(a, min_v, max_v):
        x = np.arange(min_v, max_v + 1, dtype='float')
        pmf = x ** a
        pmf /= pmf.sum()
        return stats.rv_discrete(values=(x, pmf))

    print("expected degrees: [{},{}]".format(min_degree, max_degree))
    # set in-degree and out-degree distribution
    rv_indegree = discrete_power_law(powerlaw_exponent, min_degree, max_degree)
    rv_outdegree = discrete_power_law(powerlaw_exponent, min_degree, max_degree)

    # define the return function for in and out degrees
    def degree_distribution_function(rv1, rv2):
        return (rv1.rvs(size=1), rv2.rvs(size=1))

    num_blocks, block_membership_vector = generate_block_membership(num_blocks, N_adjusted, block_size_heterogeneity)

    ####################
    # GENERATE DEGREE-CORRECTED SBM
    ####################
    blockmodel = fill_blockmodel(N_adjusted, num_blocks, ratio_within_over_between, block_membership_vector)

    if N_adjusted > 1000000:
        total_degrees = np.asarray([rv_outdegree.rvs() for i in
                                    range(N_adjusted)])
    else:
        total_degrees = rv_outdegree.rvs(size=N_adjusted)
    out_degrees = np.random.uniform(size=N_adjusted) * total_degrees
    out_degrees = np.round(out_degrees)
    in_degrees = total_degrees - out_degrees
    sum_degrees = total_degrees.sum()
    print("sum degrees: ", sum_degrees)

    expected_e = sum_degrees
    K = expected_e / (np.sum(out_degrees + in_degrees))

    print("out: [{},{}]".format(np.min(out_degrees), np.max(out_degrees)))
    print("in: [{},{}]".format(np.min(in_degrees), np.max(in_degrees)))
    # print("B:\n", blockmodel)
    graph = gt.generate_sbm(
            # Block membership of each vertex
            b=block_membership_vector,
            # Edge propensities between communities
            probs=blockmodel * (expected_e / blockmodel.sum()),
            # The out degree propensity of each vertex
            out_degs=out_degrees,
            # The in degree propensity of each vertex
            in_degs=in_degrees,
            directed=True,
            micro_ers=False,  # If True, num edges b/n groups will be exactly probs
            micro_degs=False  # If True, degrees of nodes will be exactly degs
    )

    # remove (1-density) percent of the edges
    edge_filter = graph.new_edge_property('bool')
    edge_filter.a = stats.bernoulli.rvs(density, size=edge_filter.a.shape)
    graph.set_edge_filter(edge_filter)
    graph.purge_edges()

    # remove all island vertices
    print('Filtering out zero vertices...')
    degrees = graph.get_total_degrees(np.arange(graph.num_vertices()))
    vertex_filter = graph.new_vertex_property('bool', vals=degrees > 0.0)
    graph.set_vertex_filter(vertex_filter)
    graph.purge_vertices()

    # store the nodal block memberships in a vertex property
    block_membership_vector = block_membership_vector[degrees > 0.0]
    true_partition = block_membership_vector
    assert block_membership_vector.size == graph.num_vertices()
    block_membership = graph.new_vertex_property("int", vals=block_membership_vector)

    # compute and report basic statistics on the generated graph
    analyze_generated_graph(graph, block_membership, num_blocks)

    for block_exponent in [0.25, 0.35, 0.45]:
        for block_size_variation in [1.0, 3.0, 5.0]:
            for overlap in [1.0, 3.0, 5.0]:
                # generate new blockmodel
                num_blocks = int(graph.num_vertices() ** block_exponent)
                num_blocks, block_membership_vector = generate_block_membership(num_blocks, graph.num_vertices(),
                                                                                block_size_variation)
                blockmodel = fill_blockmodel(graph.num_vertices(), num_blocks, overlap, block_membership_vector)
                block_membership = graph.new_vertex_property("int", vals=block_membership_vector)

                # shuffle the graph
                def edge_probability_fxn(row: int, col: int):
                    return blockmodel[row, col]

                rejected = gt.random_rewire(graph, "blockmodel", edge_probs=edge_probability_fxn,
                                            block_membership=block_membership)
                print("Num rejected edge moves: ", rejected)
                analyze_generated_graph(graph, block_membership, num_blocks)
                file_name = generate_filename(graph.num_vertices(), num_blocks, args.maxdegree, powerlaw_exponent,
                                              overlap, block_size_variation, density, args)
                print(file_name)
                save_graph(graph, block_membership_vector, utility, file_name)
# End of generate()


if __name__ == "__main__":
    args = parse_args()
    utility = os
    if args.remote:
        print("Getting credentials for logging into remote system")
        username = input("Username: ")
        password = getpass.getpass()
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(args.remote, username=username, password=password, timeout=2000)
        ssh.exec_command("mkdir -p {}".format(args.directory))
        utility = ssh
    generate(args, utility)
    exit()
    